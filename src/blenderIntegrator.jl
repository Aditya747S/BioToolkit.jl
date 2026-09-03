module BlenderIntegrator

using JSON
using Printf
using Statistics

export AbstractBlenderObject, AbstractBlenderPayload
export BlenderCamera, BlenderLight, BlenderMaterial, BlenderScene
export BlenderProteinPayload, BlenderSpatialPayload, BlenderContactPayload, BlenderHiCPayload, BlenderGWASPayload, BlenderRNAVelocityPayload, BlenderHMMPayload, BlenderTextPayload
export to_blender_payload, to_blender_scene, to_python_dict, write_payload_json, build_bpy_script, write_blender_script, extract_coordinates
export detect_blender_executable, is_blender_available, render_blender_headless
export export_gltf_pure, export_obj_pure
export auto_frame_camera, quick_render, open_in_blender, quick_export_3d, render_turntable

"""
    to_blender_scene(obj; kwargs...)

Generic dispatch method for creating a `BlenderScene` from domain payloads or objects.
"""
function to_blender_scene end

# -----------------------------------------------------------------------------
# 1. Scene Graph Types & Abstract Hierarchy
# -----------------------------------------------------------------------------

"""
Abstract root type for any object that can be converted into a Blender scene node.
"""
abstract type AbstractBlenderObject end

"""
Abstract root type for domain payload data structures.
"""
abstract type AbstractBlenderPayload <: AbstractBlenderObject end

"""
    BlenderCamera(; location=(0.0, -50.0, 20.0), rotation=(60.0, 0.0, 0.0), fov=45.0)

Represents camera positioning, rotation (Euler angles in degrees), and field of view.
"""
struct BlenderCamera
    location::NTuple{3, Float64}
    rotation::NTuple{3, Float64} # Euler degrees
    fov::Float64
end
BlenderCamera(; location=(0.0, -50.0, 20.0), rotation=(60.0, 0.0, 0.0), fov=45.0) = 
    BlenderCamera(location, rotation, fov)

"""
    BlenderLight(type::Symbol, location::NTuple{3,Float64}, energy::Real; color=(1.0, 1.0, 1.0))

Represents a light source (:area, :sun, :point, :spot).
"""
struct BlenderLight
    light_type::Symbol
    location::NTuple{3, Float64}
    energy::Float64
    color::NTuple{3, Float64}
end
BlenderLight(type::Symbol, loc::NTuple{3, Float64}, energy::Real; color=(1.0, 1.0, 1.0)) = 
    BlenderLight(type, loc, Float64(energy), color)

"""
    BlenderMaterial(; name="BioMaterial", color=(0.8, 0.8, 0.8, 1.0), metallic=0.0, roughness=0.4, transmission=0.0, emission_energy=0.0)

Material properties for PBR shader rendering (Cycles / EEVEE).
"""
struct BlenderMaterial
    name::String
    color::NTuple{4, Float64} # RGBA
    metallic::Float64
    roughness::Float64
    transmission::Float64 # glassmorphism
    emission_energy::Float64
end
function BlenderMaterial(; 
    name::String="BioMaterial",
    color::NTuple{4, Float64}=(0.8, 0.8, 0.8, 1.0),
    metallic::Real=0.0,
    roughness::Real=0.4,
    transmission::Real=0.0,
    emission_energy::Real=0.0
)
    return BlenderMaterial(name, color, Float64(metallic), Float64(roughness), Float64(transmission), Float64(emission_energy))
end

"""
    BlenderScene(; objects=AbstractBlenderObject[], camera=BlenderCamera(), lights=[...], engine=:cycles, samples=128, resolution=(1920, 1080))

Top-level Blender scene container.
"""
struct BlenderScene
    objects::Vector{AbstractBlenderObject}
    camera::BlenderCamera
    lights::Vector{BlenderLight}
    engine::Symbol # :cycles, :eevee
    samples::Int
    resolution::Tuple{Int, Int}
end

function BlenderScene(; 
    objects=AbstractBlenderObject[], 
    camera=BlenderCamera(), 
    lights=[BlenderLight(:area, (10.0, -10.0, 30.0), 1000.0)],
    engine=:cycles,
    samples=128,
    resolution=(1920, 1080)
)
    return BlenderScene(objects, camera, lights, engine, samples, resolution)
end

function to_blender_scene(obj::AbstractBlenderObject; camera=nothing, lights=nothing, engine=:cycles, samples=128, resolution=(1920, 1080))
    cam = camera === nothing ? auto_frame_camera(obj) : camera
    lts = lights === nothing ? [BlenderLight(:area, (10.0, -10.0, 30.0), 1000.0)] : lights
    return BlenderScene(objects=[obj], camera=cam, lights=lts, engine=engine, samples=samples, resolution=resolution)
end

function to_blender_scene(domain_obj; kwargs...)
    payload = to_blender_payload(domain_obj)
    return to_blender_scene(payload; kwargs...)
end

# -----------------------------------------------------------------------------
# 2. Generic Extension & Serialization Hooks
# -----------------------------------------------------------------------------

"""
    to_blender_payload(obj; kwargs...)

Generic dispatch method. Extended by BioToolkit domain modules (structure, spatial, hic, etc.).
"""
function to_blender_payload end

"""
    to_python_dict(obj::AbstractBlenderObject) -> Dict{String, Any}

Converts a Blender payload struct into a dictionary ready for JSON bootstrap export.
"""
function to_python_dict end

function material_to_dict(mat::BlenderMaterial)
    return Dict{String, Any}(
        "name" => mat.name,
        "color" => collect(mat.color),
        "metallic" => mat.metallic,
        "roughness" => mat.roughness,
        "transmission" => mat.transmission,
        "emission_energy" => mat.emission_energy
    )
end

# -----------------------------------------------------------------------------
# 3. Domain Payload Structs & Serialization
# -----------------------------------------------------------------------------

"""
    BlenderProteinPayload

Atomic coordinates and secondary structure styling.
"""
struct BlenderProteinPayload <: AbstractBlenderPayload
    name::String
    coordinates::Matrix{Float64} # N x 3
    elements::Vector{String}
    residue_indices::Vector{Int}
    style::Symbol # :cartoon, :vdw, :ball_and_stick
    color_scheme::Symbol # :element, :rainbow, :sasa, :secondary_structure
    material::BlenderMaterial
    scores::Vector{Float64}
end

function BlenderProteinPayload(name::String, coordinates::Matrix{Float64}, elements::Vector{String}, residue_indices::Vector{Int}, style::Symbol, color_scheme::Symbol, material::BlenderMaterial; scores::Vector{Float64}=Float64[])
    return BlenderProteinPayload(name, coordinates, elements, residue_indices, style, color_scheme, material, scores)
end

to_python_dict(dict::AbstractDict) = Dict{String, Any}(string(k) => v for (k, v) in dict)

function to_python_dict(obj::BlenderProteinPayload)
    return Dict{String, Any}(
        "type" => "protein",
        "geometry" => "protein",
        "name" => obj.name,
        "coordinates" => [vec(obj.coordinates[i, :]) for i in 1:size(obj.coordinates, 1)],
        "elements" => obj.elements,
        "residue_indices" => obj.residue_indices,
        "style" => string(obj.style),
        "color_scheme" => string(obj.color_scheme),
        "scores" => obj.scores,
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderSpatialPayload

3D spatial transcriptomics spot coordinates and cluster colors.
"""
struct BlenderSpatialPayload <: AbstractBlenderPayload
    name::String
    coordinates::Matrix{Float64} # N x 3
    cluster_labels::Vector{String}
    colors::Vector{NTuple{3, Float64}}
    glyph_scale::Float64
    material::BlenderMaterial
end

function to_python_dict(obj::BlenderSpatialPayload)
    return Dict{String, Any}(
        "type" => "spatial",
        "geometry" => "points",
        "name" => obj.name,
        "coordinates" => [vec(obj.coordinates[i, :]) for i in 1:size(obj.coordinates, 1)],
        "cluster_labels" => obj.cluster_labels,
        "colors" => [collect(c) for c in obj.colors],
        "glyph_scale" => obj.glyph_scale,
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderContactPayload

Coevolutionary Direct Information (DCA) contact pairs and energy tube connections.
"""
struct BlenderContactPayload <: AbstractBlenderPayload
    name::String
    coordinates::Matrix{Float64} # N x 3
    contact_pairs::Vector{Tuple{Int, Int}}
    coupling_scores::Vector{Float64}
    tube_radius::Float64
    material::BlenderMaterial
end

function to_python_dict(obj::BlenderContactPayload)
    return Dict{String, Any}(
        "type" => "contact_network",
        "geometry" => "curve",
        "name" => obj.name,
        "coordinates" => [vec(obj.coordinates[i, :]) for i in 1:size(obj.coordinates, 1)],
        "contact_pairs" => [[p[1], p[2]] for p in obj.contact_pairs],
        "coupling_scores" => obj.coupling_scores,
        "tube_radius" => obj.tube_radius,
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderHiCPayload

3D reconstructed chromatin folding splines with Bevel tubes and TAD boundary anchors.
"""
struct BlenderHiCPayload <: AbstractBlenderPayload
    name::String
    spline_points::Matrix{Float64} # N x 3
    compartment_scores::Vector{Float64}
    loop_anchors::Vector{Tuple{Int, Int}}
    tube_radius::Float64
    material::BlenderMaterial
end

function to_python_dict(obj::BlenderHiCPayload)
    return Dict{String, Any}(
        "type" => "hic_spline",
        "geometry" => "curve",
        "name" => obj.name,
        "spline_points" => [vec(obj.spline_points[i, :]) for i in 1:size(obj.spline_points, 1)],
        "compartment_scores" => obj.compartment_scores,
        "loop_anchors" => [[a[1], a[2]] for a in obj.loop_anchors],
        "tube_radius" => obj.tube_radius,
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderGWASPayload

3D topographic terrain map (X=chrom, Y=trait, Z=-log10(p)).
"""
struct BlenderGWASPayload <: AbstractBlenderPayload
    name::String
    height_matrix::Matrix{Float64} # M x N heightmap
    chromosome_labels::Vector{String}
    trait_labels::Vector{String}
    material::BlenderMaterial
end

function to_python_dict(obj::BlenderGWASPayload)
    return Dict{String, Any}(
        "type" => "gwas_terrain",
        "geometry" => "mesh",
        "name" => obj.name,
        "height_matrix" => [vec(obj.height_matrix[i, :]) for i in 1:size(obj.height_matrix, 1)],
        "chromosome_labels" => obj.chromosome_labels,
        "trait_labels" => obj.trait_labels,
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderHMMPayload

3D HMM state transition topology network / posterior probability landscape for Blender 3D.
"""
struct BlenderHMMPayload <: AbstractBlenderPayload
    name::String
    states::Vector{String}
    transition_matrix::Matrix{Float64} # N x N probability matrix
    posterior_matrix::Matrix{Float64}   # N x L state probability heatmap matrix
    material::BlenderMaterial
end

function to_python_dict(obj::BlenderHMMPayload)
    return Dict{String, Any}(
        "type" => "hmm_landscape",
        "geometry" => "hmm_landscape",
        "name" => obj.name,
        "states" => obj.states,
        "transition_matrix" => [vec(obj.transition_matrix[i, :]) for i in 1:size(obj.transition_matrix, 1)],
        "posterior_matrix" => [vec(obj.posterior_matrix[i, :]) for i in 1:size(obj.posterior_matrix, 1)],
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderRNAVelocityPayload

3D vector field representing cell differentiation state flows.
"""
struct BlenderRNAVelocityPayload <: AbstractBlenderPayload
    name::String
    origins::Matrix{Float64} # N x 3
    vectors::Matrix{Float64} # N x 3
    arrow_scale::Float64
    material::BlenderMaterial
end

function to_python_dict(obj::BlenderRNAVelocityPayload)
    return Dict{String, Any}(
        "type" => "rna_velocity",
        "geometry" => "vector_field",
        "name" => obj.name,
        "origins" => [vec(obj.origins[i, :]) for i in 1:size(obj.origins, 1)],
        "vectors" => [vec(obj.vectors[i, :]) for i in 1:size(obj.vectors, 1)],
        "arrow_scale" => obj.arrow_scale,
        "material" => material_to_dict(obj.material)
    )
end

"""
    BlenderTextPayload(text::String; name="TextLabel", location=(0.0, 0.0, 0.0), scale=1.0, color=(1.0, 1.0, 1.0, 1.0))

Represents 3D text annotation labels inside the Blender scene graph.
"""
struct BlenderTextPayload <: AbstractBlenderPayload
    name::String
    text::String
    location::NTuple{3, Float64}
    scale::Float64
    material::BlenderMaterial
end

function BlenderTextPayload(text::String; name::String="TextLabel", location::NTuple{3, Float64}=(0.0, 0.0, 0.0), scale::Real=1.0, color::NTuple{4, Float64}=(1.0, 1.0, 1.0, 1.0))
    mat = BlenderMaterial(name=name * "_mat", color=color)
    return BlenderTextPayload(name, text, location, Float64(scale), mat)
end

function to_python_dict(obj::BlenderTextPayload)
    return Dict{String, Any}(
        "type" => "text",
        "geometry" => "text",
        "name" => obj.name,
        "text" => obj.text,
        "location" => collect(obj.location),
        "scale" => obj.scale,
        "material" => material_to_dict(obj.material)
    )
end

# -----------------------------------------------------------------------------
# 4. System Binary Autodetection
# -----------------------------------------------------------------------------

"""
    detect_blender_executable()

Locates the `blender` binary on the system PATH or via `ENV["BIOTOOLKIT_BLENDER_PATH"]`.
Returns `nothing` if Blender is not installed (non-blocking).
"""
function detect_blender_executable()
    env_path = get(ENV, "BIOTOOLKIT_BLENDER_PATH", "")
    if !isempty(env_path) && isfile(env_path)
        return env_path
    end
    return Sys.which("blender")
end

is_blender_available() = detect_blender_executable() !== nothing

# -----------------------------------------------------------------------------
# 5. Bootstrap + JSON Payload Serialization
# -----------------------------------------------------------------------------

"""
    write_payload_json(scene::BlenderScene, path::String)

Serializes all objects in `scene` to a clean JSON file for Blender python bootstrap loading.
Uses standard JSON library.
"""
function write_payload_json(scene::BlenderScene, path::String)
    abs_path = abspath(path)
    payload_dict = Dict{String, Any}(
        "camera" => Dict(
            "location" => collect(scene.camera.location),
            "rotation" => collect(scene.camera.rotation),
            "fov" => scene.camera.fov
        ),
        "lights" => [
            Dict(
                "type" => string(l.light_type),
                "location" => collect(l.location),
                "energy" => l.energy,
                "color" => collect(l.color)
            ) for l in scene.lights
        ],
        "engine" => string(scene.engine),
        "samples" => scene.samples,
        "resolution" => collect(scene.resolution),
        "objects" => [to_python_dict(obj) for obj in scene.objects]
    )
    
    json_str = JSON.json(payload_dict)
    write(abs_path, json_str)
    return abs_path
end

# -----------------------------------------------------------------------------
# 6. Static Python Bootstrap Script Generator
# -----------------------------------------------------------------------------

"""
    build_bpy_script(scene::BlenderScene; data_path="", output_path="", format=:png) -> String

Formats a clean Python script string that loads `data_path` JSON payload and renders/exports.
"""
function build_bpy_script(scene::BlenderScene; data_path::String="", output_path::String="", format::Symbol=:png, turntable::Bool=false, frames::Int=36, fps::Int=24)
    abs_data_path = isempty(data_path) ? "" : abspath(data_path)
    abs_output_path = isempty(output_path) ? "" : abspath(output_path)

    json_data_path_literal = JSON.json(abs_data_path)
    json_output_path_literal = JSON.json(abs_output_path)
    format_str = lowercase(string(format))

    cam = scene.camera
    cam_loc = [cam.location[1], cam.location[2], cam.location[3]]
    cam_rot = [cam.rotation[1], cam.rotation[2], cam.rotation[3]]
    cam_fov = cam.fov

    engine_str = uppercase(string(scene.engine))
    samples_val = scene.samples
    res_x = scene.resolution[1]
    res_y = scene.resolution[2]

    return """
# Static Python Bootstrap Script generated by BioToolkit.jl BlenderIntegrator
import bpy
import math
import json
import os
import colorsys

def apply_material(obj, mat_data):
    if not obj or not mat_data:
        return None
    mat_name = mat_data.get('name', 'BioMat')
    mat = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    bsdf = nodes.get('Principled BSDF')
    
    c = mat_data.get('color', [0.8, 0.8, 0.8, 1.0])
    
    if hasattr(obj.data, 'color_attributes') and len(obj.data.color_attributes) > 0:
        try:
            attr_node = nodes.new(type='ShaderNodeAttribute')
            attr_node.attribute_name = obj.data.color_attributes[0].name
            if bsdf:
                links.new(attr_node.outputs['Color'], bsdf.inputs['Base Color'])
        except Exception:
            pass
    elif bsdf:
        if 'Base Color' in bsdf.inputs:
            bsdf.inputs['Base Color'].default_value = (c[0], c[1], c[2], c[3])
        if 'Metallic' in bsdf.inputs:
            bsdf.inputs['Metallic'].default_value = mat_data.get('metallic', 0.0)
        if 'Roughness' in bsdf.inputs:
            bsdf.inputs['Roughness'].default_value = mat_data.get('roughness', 0.4)
            
        trans_val = mat_data.get('transmission', 0.0)
        if 'Transmission Weight' in bsdf.inputs:
            bsdf.inputs['Transmission Weight'].default_value = trans_val
        elif 'Transmission' in bsdf.inputs:
            bsdf.inputs['Transmission'].default_value = trans_val
            
        emit_val = mat_data.get('emission_energy', 0.0)
        if 'Emission Color' in bsdf.inputs:
            bsdf.inputs['Emission Color'].default_value = (c[0], c[1], c[2])
        elif 'Emission' in bsdf.inputs:
            bsdf.inputs['Emission'].default_value = (c[0], c[1], c[2])
        if 'Emission Strength' in bsdf.inputs:
            bsdf.inputs['Emission Strength'].default_value = emit_val
            
    if hasattr(obj.data, 'materials'):
        if obj.data.materials:
            obj.data.materials[0] = mat
        else:
            obj.data.materials.append(mat)
    return mat

# Reset default factory scene
bpy.ops.wm.read_factory_settings(use_empty=True)

# Data file loading
data_file = $json_data_path_literal
payload = {}
if data_file and os.path.exists(data_file):
    with open(data_file, 'r') as f:
        payload = json.load(f)

# Engine configuration
engine_name = payload.get('engine', '$engine_str').upper()
if engine_name in ('CYCLES', 'CYCLES_GPU'):
    bpy.context.scene.render.engine = 'CYCLES'
    bpy.context.scene.cycles.samples = payload.get('samples', $samples_val)
    try:
        cycles_prefs = bpy.context.preferences.addons['cycles'].preferences
        cycles_prefs.get_devices()
        for device in cycles_prefs.devices:
            device.use = True
        bpy.context.scene.cycles.device = 'GPU'
    except Exception:
        pass
else:
    try:
        bpy.context.scene.render.engine = 'BLENDER_EEVEE_NEXT'
    except Exception:
        try:
            bpy.context.scene.render.engine = 'BLENDER_EEVEE'
        except Exception:
            pass

res = payload.get('resolution', [$res_x, $res_y])
bpy.context.scene.render.resolution_x = res[0]
bpy.context.scene.render.resolution_y = res[1]

# Camera setup
cam_cfg = payload.get('camera', {})
cam_loc = cam_cfg.get('location', $cam_loc)
cam_rot = cam_cfg.get('rotation', $cam_rot)
cam_fov = cam_cfg.get('fov', $cam_fov)

cam_data = bpy.data.cameras.new(name='BioCamera')
cam_obj = bpy.data.objects.new('BioCamera', cam_data)
bpy.context.scene.collection.objects.link(cam_obj)
bpy.context.scene.camera = cam_obj
cam_obj.location = (cam_loc[0], cam_loc[1], cam_loc[2])
cam_obj.rotation_euler = (math.radians(cam_rot[0]), math.radians(cam_rot[1]), math.radians(cam_rot[2]))
try:
    cam_data.lens_unit = 'FOV'
    cam_data.angle = math.radians(cam_fov)
except Exception:
    pass

# Lights setup
lights_data = payload.get('lights', [])
if not lights_data:
    lights_data = [{'type': 'AREA', 'location': [10.0, -10.0, 30.0], 'energy': 1000.0, 'color': [1.0, 1.0, 1.0]}]

for idx, l in enumerate(lights_data):
    l_type = str(l.get('type', 'AREA')).upper()
    l_loc = l.get('location', [10.0, -10.0, 30.0])
    l_energy = float(l.get('energy', 1000.0))
    l_color = l.get('color', [1.0, 1.0, 1.0])
    
    light_data = bpy.data.lights.new(name=f'BioLight_{idx}', type=l_type)
    light_data.energy = l_energy
    if hasattr(light_data, 'color') and len(l_color) >= 3:
        light_data.color = (float(l_color[0]), float(l_color[1]), float(l_color[2]))
    light_obj = bpy.data.objects.new(f'BioLight_{idx}', light_data)
    bpy.context.scene.collection.objects.link(light_obj)
    light_obj.location = (float(l_loc[0]), float(l_loc[1]), float(l_loc[2]))

# Object payloads processing
for idx, obj in enumerate(payload.get('objects', [])):
    obj_type = obj.get('type', '')
    geom = obj.get('geometry', '')
    if not geom:
        if obj_type in ('protein',): geom = 'protein'
        elif obj_type in ('spatial',): geom = 'points'
        elif obj_type in ('contact_network', 'hic_spline'): geom = 'curve'
        elif obj_type in ('gwas_terrain', 'hmm_landscape'): geom = 'mesh'
        elif obj_type in ('rna_velocity',): geom = 'vector_field'
        elif obj_type in ('text',): geom = 'text'
        else: geom = 'points'
        
    name = obj.get('name', f'BioObj_{idx}')
    coords = obj.get('coordinates', [])
    mat_data = obj.get('material')
    created_obj = None

    if geom in ('points', 'spatial', 'protein'):
        mesh = bpy.data.meshes.new(name=f'{name}_mesh')
        created_obj = bpy.data.objects.new(name, mesh)
        bpy.context.scene.collection.objects.link(created_obj)
        mesh.from_pydata(coords, [], [])
        mesh.update()
        
        elements = obj.get('elements', [])
        color_scheme = obj.get('color_scheme', 'element')
        scores = obj.get('scores', [])
        n_verts = len(mesh.vertices)
        colors = []
        
        if color_scheme == 'element' and elements:
            cpk_map = {'C': (0.5, 0.5, 0.5, 1.0), 'N': (0.0, 0.0, 1.0, 1.0), 'O': (1.0, 0.0, 0.0, 1.0), 'S': (1.0, 1.0, 0.0, 1.0), 'H': (1.0, 1.0, 1.0, 1.0)}
            colors = [cpk_map.get(str(el).upper(), (0.8, 0.8, 0.8, 1.0)) for el in elements]
        elif color_scheme == 'rainbow':
            for i in range(n_verts):
                frac = i / max(n_verts - 1, 1)
                hue = (1.0 - frac) * 0.666
                rgb = colorsys.hsv_to_rgb(hue, 0.9, 0.95)
                colors.append((rgb[0], rgb[1], rgb[2], 1.0))
        elif color_scheme == 'secondary_structure':
            res_indices = obj.get('residue_indices', [])
            for i in range(n_verts):
                res_idx = res_indices[i] if i < len(res_indices) else i
                mod = (res_idx // 10) % 3
                if mod == 0:
                    colors.append((0.9, 0.2, 0.2, 1.0))
                elif mod == 1:
                    colors.append((0.9, 0.8, 0.1, 1.0))
                else:
                    colors.append((0.2, 0.8, 0.4, 1.0))
        elif color_scheme == 'sasa' or (color_scheme in ('b_factor', 'score') and scores):
            min_s, max_s = (min(scores), max(scores)) if scores else (0.0, 1.0)
            rng_s = max(max_s - min_s, 1e-5)
            for i in range(n_verts):
                s_val = scores[i] if i < len(scores) else min_s
                norm_s = (s_val - min_s) / rng_s
                if color_scheme == 'sasa':
                    r_col = min(1.0, 2.0 * norm_s)
                    g_col = min(1.0, 2.0 * (1.0 - abs(norm_s - 0.5)))
                    b_col = min(1.0, 2.0 * (1.0 - norm_s))
                    colors.append((r_col, g_col, b_col, 1.0))
                else:
                    hue = (1.0 - norm_s) * 0.666
                    rgb = colorsys.hsv_to_rgb(hue, 0.95, 0.95)
                    colors.append((rgb[0], rgb[1], rgb[2], 1.0))
                    
        if len(colors) == n_verts:
            try:
                color_attr = mesh.color_attributes.new(name='Color', type='FLOAT_COLOR', domain='POINT')
                for i, col in enumerate(colors):
                    color_attr.data[i].color = col
            except Exception:
                pass
                
        # Sphere instancing for vertices so points/atoms render 3D geometry
        glyph_scale = float(obj.get('glyph_scale', 0.5))
        if geom == 'protein':
            style = obj.get('style', 'vdw')
            if style == 'cartoon':
                glyph_scale = 0.3
            elif style == 'ball_and_stick':
                glyph_scale = 0.25
            elif style == 'vdw':
                glyph_scale = 0.8
                
        try:
            bpy.ops.mesh.primitive_uv_sphere_add(radius=glyph_scale, location=(0, 0, 0))
            sph_obj = bpy.context.active_object
            sph_obj.name = f'{name}_sphere'
            sph_obj.parent = created_obj
            created_obj.instance_type = 'VERTS'
            created_obj.show_instancer_for_render = False
            if mat_data:
                apply_material(sph_obj, mat_data)
        except Exception:
            pass
            
        # Additional backbone trace for protein cartoon/ball_and_stick
        if geom == 'protein' and obj.get('style') in ('cartoon', 'ball_and_stick') and len(coords) > 1:
            try:
                bb_curve = bpy.data.curves.new(f'{name}_trace', type='CURVE')
                bb_curve.dimensions = '3D'
                bb_curve.bevel_depth = 0.15 if obj.get('style') == 'cartoon' else 0.08
                bb_obj = bpy.data.objects.new(f'{name}_trace', bb_curve)
                bpy.context.scene.collection.objects.link(bb_obj)
                spline = bb_curve.splines.new('POLY')
                spline.points.add(len(coords) - 1)
                for i_c, p_c in enumerate(coords):
                    spline.points[i_c].co = (p_c[0], p_c[1], p_c[2], 1)
                if mat_data:
                    apply_material(bb_obj, mat_data)
            except Exception:
                pass

    elif geom == 'curve':
        curve = bpy.data.curves.new(f'{name}_curve', type='CURVE')
        curve.dimensions = '3D'
        curve.bevel_depth = float(obj.get('tube_radius', 0.2))
        created_obj = bpy.data.objects.new(name, curve)
        bpy.context.scene.collection.objects.link(created_obj)
        if obj_type == 'contact_network':
            pairs = obj.get('contact_pairs', [])
            coupling = obj.get('coupling_scores', [])
            min_c = min(coupling) if coupling else 0.0
            max_c = max(coupling) if coupling else 1.0
            rng_c = max(max_c - min_c, 1e-5)
            for p_idx, p in enumerate(pairs):
                spline = curve.splines.new('POLY')
                spline.points.add(1)
                p1, p2 = coords[p[0]-1], coords[p[1]-1]
                spline.points[0].co = (p1[0], p1[1], p1[2], 1)
                spline.points[1].co = (p2[0], p2[1], p2[2], 1)
                if coupling and p_idx < len(coupling):
                    c_norm = (coupling[p_idx] - min_c) / rng_c
                    r_val = 0.5 + 1.5 * c_norm
                    spline.points[0].radius = r_val
                    spline.points[1].radius = r_val
        else:
            pts = obj.get('spline_points', [])
            comp_scores = obj.get('compartment_scores', [])
            if len(pts) > 0:
                spline = curve.splines.new('NURBS')
                spline.use_endpoint_u = True
                spline.points.add(len(pts)-1)
                min_cs = min(comp_scores) if comp_scores else -1.0
                max_cs = max(comp_scores) if comp_scores else 1.0
                rng_cs = max(max_cs - min_cs, 1e-5)
                for i, p in enumerate(pts):
                    spline.points[i].co = (p[0], p[1], p[2], 1)
                    if comp_scores and i < len(comp_scores):
                        cs_norm = (comp_scores[i] - min_cs) / rng_cs
                        spline.points[i].radius = 0.5 + 1.0 * cs_norm

    elif geom in ('mesh', 'hmm_landscape'):
        h_matrix = obj.get('height_matrix', [])
        trans_matrix = obj.get('transition_matrix', [])
        post_matrix = obj.get('posterior_matrix', [])
        
        if not h_matrix and post_matrix:
            h_matrix = [[float(val) * 10.0 for val in row] for row in post_matrix]

        rows = len(h_matrix)
        cols = len(h_matrix[0]) if rows > 0 else 0
        verts = obj.get('vertices', [])
        faces = obj.get('faces', [])
        if h_matrix and not verts and rows > 0 and cols > 0:
            verts = [(float(x), float(y), float(h_matrix[x][y])) for x in range(rows) for y in range(cols)]
            faces = []
            for x in range(rows - 1):
                for y in range(cols - 1):
                    i1 = x * cols + y
                    i2 = x * cols + y + 1
                    i3 = (x + 1) * cols + y + 1
                    i4 = (x + 1) * cols + y
                    faces.append([i1, i2, i3, i4])
        if verts:
            mesh = bpy.data.meshes.new(f'{name}_mesh')
            mesh.from_pydata(verts, [], faces)
            mesh.update()
            if h_matrix:
                heights = [v[2] for v in verts]
                min_h, max_h = min(heights), max(heights)
                rng_h = max(max_h - min_h, 1e-5)
                try:
                    color_attr = mesh.color_attributes.new(name='Color', type='FLOAT_COLOR', domain='POINT')
                    for i, h in enumerate(heights):
                        norm_h = (h - min_h) / rng_h
                        r = min(1.0, max(0.0, 1.5 * norm_h - 0.2))
                        g = min(1.0, max(0.0, 2.0 * norm_h - 1.0))
                        b = min(1.0, max(0.0, 1.0 - 1.5 * norm_h))
                        color_attr.data[i].color = (r, g, b, 1.0)
                except Exception:
                    pass
            created_obj = bpy.data.objects.new(name, mesh)
            bpy.context.scene.collection.objects.link(created_obj)
            
        if trans_matrix:
            n_states = len(trans_matrix)
            states = obj.get('states', [f'S{i}' for i in range(n_states)])
            state_coords = []
            radius = max(float(n_states) * 1.5, 5.0)
            for i in range(n_states):
                angle = 2.0 * math.pi * i / max(n_states, 1)
                sx = radius * math.cos(angle)
                sy = radius * math.sin(angle)
                sz = 0.0
                state_coords.append((sx, sy, sz))
                try:
                    bpy.ops.mesh.primitive_uv_sphere_add(radius=0.8, location=(sx, sy, sz))
                    st_obj = bpy.context.active_object
                    st_obj.name = f'{name}_State_{states[i]}'
                    if mat_data:
                        apply_material(st_obj, mat_data)
                except Exception:
                    pass
                    
            tr_curve = bpy.data.curves.new(f'{name}_trans_curve', type='CURVE')
            tr_curve.dimensions = '3D'
            tr_curve.bevel_depth = 0.15
            tr_obj = bpy.data.objects.new(f'{name}_transitions', tr_curve)
            bpy.context.scene.collection.objects.link(tr_obj)
            if not created_obj:
                created_obj = tr_obj
            for i in range(n_states):
                for j in range(n_states):
                    prob = float(trans_matrix[i][j])
                    if prob > 0.01:
                        spline = tr_curve.splines.new('POLY')
                        spline.points.add(1)
                        p1, p2 = state_coords[i], state_coords[j]
                        spline.points[0].co = (p1[0], p1[1], p1[2], 1)
                        spline.points[1].co = (p2[0], p2[1], p2[2], 1)
                        spline.points[0].radius = 0.2 + 0.8 * prob
                        spline.points[1].radius = 0.2 + 0.8 * prob

    elif geom == 'vector_field':
        origins = obj.get('origins', [])
        vectors = obj.get('vectors', [])
        arrow_scale = float(obj.get('arrow_scale', 1.0))
        curve = bpy.data.curves.new(f'{name}_velo', type='CURVE')
        curve.dimensions = '3D'
        curve.bevel_depth = 0.05 * arrow_scale
        created_obj = bpy.data.objects.new(name, curve)
        bpy.context.scene.collection.objects.link(created_obj)
        for o, v in zip(origins, vectors):
            spline = curve.splines.new('POLY')
            spline.points.add(1)
            spline.points[0].co = (o[0], o[1], o[2], 1)
            spline.points[1].co = (o[0] + v[0]*arrow_scale, o[1] + v[1]*arrow_scale, o[2] + v[2]*arrow_scale, 1)
            
            # Arrow cone head
            tip_x = o[0] + v[0]*arrow_scale
            tip_y = o[1] + v[1]*arrow_scale
            tip_z = o[2] + v[2]*arrow_scale
            try:
                bpy.ops.mesh.primitive_cone_add(radius1=0.15*arrow_scale, depth=0.3*arrow_scale, location=(tip_x, tip_y, tip_z))
                c_obj = bpy.context.active_object
                if mat_data:
                    apply_material(c_obj, mat_data)
            except Exception:
                pass

    elif geom == 'text':
        text_str = obj.get('text', '')
        font_curve = bpy.data.curves.new(name=f'{name}_font', type='FONT')
        font_curve.body = text_str
        created_obj = bpy.data.objects.new(name, font_curve)
        bpy.context.scene.collection.objects.link(created_obj)
        loc = obj.get('location', [0.0, 0.0, 0.0])
        created_obj.location = (loc[0], loc[1], loc[2])
        scale = obj.get('scale', 1.0)
        created_obj.scale = (scale, scale, scale)

    if created_obj and mat_data:
        apply_material(created_obj, mat_data)

# Turntable keyframe animation
if $(turntable ? "True" : "False"):
    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = $frames
    scene.render.fps = $fps
    
    # Calculate scene centroid dynamically
    all_objs = [o for o in bpy.context.scene.collection.objects if o.type in ('MESH', 'CURVE', 'FONT')]
    if all_objs:
        target_x = sum(o.location.x for o in all_objs) / len(all_objs)
        target_y = sum(o.location.y for o in all_objs) / len(all_objs)
        target_z = sum(o.location.z for o in all_objs) / len(all_objs)
    else:
        target_x, target_y, target_z = 0.0, 0.0, 0.0

    dx_init = cam_obj.location.x - target_x
    dy_init = cam_obj.location.y - target_y
    dist_r = math.sqrt(dx_init**2 + dy_init**2)
    if dist_r < 1e-3:
        dist_r = 40.0
    cam_z = cam_obj.location.z

    for f in range(1, $frames + 1):
        scene.frame_set(f)
        angle = 2.0 * math.pi * (f - 1) / $frames
        cam_obj.location.x = target_x + dist_r * math.sin(angle)
        cam_obj.location.y = target_y - dist_r * math.cos(angle)
        cam_obj.location.z = cam_z
        
        dx = target_x - cam_obj.location.x
        dy = target_y - cam_obj.location.y
        dz = target_z - cam_obj.location.z
        yaw = math.atan2(dx, -dy)
        pitch = math.atan2(dz, math.sqrt(dx*dx + dy*dy))
        cam_obj.rotation_euler = (math.pi/2 - pitch, 0, yaw)
        cam_obj.keyframe_insert(data_path='location', index=-1)
        cam_obj.keyframe_insert(data_path='rotation_euler', index=-1)

# Render & export output handler
if $json_output_path_literal:
    out_file = $json_output_path_literal
    fmt = "$format_str"
    scene = bpy.context.scene
    scene.render.filepath = out_file

    if fmt in ('mp4', 'avi'):
        scene.render.image_settings.file_format = 'FFMPEG'
        scene.render.ffmpeg.format = 'MPEG4'
        scene.render.ffmpeg.codec = 'H264'
        if not $(turntable ? "True" : "False"):
            scene.frame_start = 1
            scene.frame_end = 1
        bpy.ops.render.render(animation=True)
    elif fmt == 'gif':
        scene.render.image_settings.file_format = 'FFMPEG'
        try:
            scene.render.ffmpeg.format = 'GIF'
        except Exception:
            scene.render.ffmpeg.format = 'MPEG4'
        if not $(turntable ? "True" : "False"):
            scene.frame_start = 1
            scene.frame_end = 1
        bpy.ops.render.render(animation=True)
    elif fmt in ('jpg', 'jpeg'):
        scene.render.image_settings.file_format = 'JPEG'
        if $(turntable ? "True" : "False"):
            bpy.ops.render.render(animation=True)
        else:
            bpy.ops.render.render(write_still=True)
    elif fmt == 'png':
        scene.render.image_settings.file_format = 'PNG'
        if $(turntable ? "True" : "False"):
            bpy.ops.render.render(animation=True)
        else:
            bpy.ops.render.render(write_still=True)
    elif fmt in ('gltf', 'glb'):
        if hasattr(bpy.ops.export_scene, 'gltf'):
            bpy.ops.export_scene.gltf(filepath=out_file, export_format='GLB')
        elif hasattr(bpy.ops.wm, 'gltf_export'):
            bpy.ops.wm.gltf_export(filepath=out_file)
    elif fmt in ('usd', 'usda', 'usdc'):
        bpy.ops.wm.usd_export(filepath=out_file)
    elif fmt == 'obj':
        if hasattr(bpy.ops.wm, 'obj_export'):
            bpy.ops.wm.obj_export(filepath=out_file)
        elif hasattr(bpy.ops.export_scene, 'obj'):
            bpy.ops.export_scene.obj(filepath=out_file)
"""
end

"""
    write_blender_script(scene::BlenderScene, path::String; data_path="")

Writes the generated Blender Python script to disk.
"""
function write_blender_script(scene::BlenderScene, path::String; data_path::String="")
    abs_path = abspath(path)
    script_str = build_bpy_script(scene; data_path=data_path)
    write(abs_path, script_str)
    return abs_path
end

"""
    render_blender_headless(scene::BlenderScene, output_path::String; format::Symbol=:png, timeout::Int=300, turntable::Bool=false, frames::Int=36, fps::Int=24)

Launches Blender in background mode using standard Julia OS process pipes `run()`.
Safely handles process execution with output verification and stderr capture.
"""
function render_blender_headless(scene::BlenderScene, output_path::String; format::Symbol=:png, timeout::Int=300, turntable::Bool=false, frames::Int=36, fps::Int=24)
    blender_bin = detect_blender_executable()
    if blender_bin === nothing
        error("Blender executable not found on system PATH. Install Blender or set ENV[\"BIOTOOLKIT_BLENDER_PATH\"].")
    end

    abs_output = abspath(output_path)
    temp_json = tempname() * ".json"
    write_payload_json(scene, temp_json)

    script_str = build_bpy_script(scene; data_path=temp_json, output_path=abs_output, format=format, turntable=turntable, frames=frames, fps=fps)
    temp_script = tempname() * ".py"
    write(temp_script, script_str)

    effective_timeout = turntable ? max(timeout, frames * 10) : timeout

    out_buf = IOBuffer()
    err_buf = IOBuffer()

    try
        cmd = `$blender_bin --background --python $temp_script`
        proc = run(pipeline(cmd, stdout=out_buf, stderr=err_buf), wait=false)
        t0 = time()
        while process_running(proc)
            if time() - t0 > effective_timeout
                kill(proc)
                wait(proc)
                err_text = String(take!(err_buf))
                error("Blender execution timed out after $(effective_timeout)s.\nStderr:\n$err_text")
            end
            sleep(0.1)
        end
        wait(proc)

        err_text = String(take!(err_buf))
        out_text = String(take!(out_buf))

        if proc.exitcode != 0
            error("Blender process exited with error code $(proc.exitcode).\nStderr:\n$err_text\nStdout:\n$out_text")
        end

        # Verify rendered output exists and is non-empty
        if !isfile(abs_output) || filesize(abs_output) == 0
            # Handle turntable image sequence naming if applicable (e.g. frame 0001)
            if turntable && (format in (:png, :jpg, :jpeg))
                base, ext = splitext(abs_output)
                frame1_path = "$(base)0001$(ext)"
                if isfile(frame1_path) && filesize(frame1_path) > 0
                    return frame1_path
                end
            end
            error("Blender executed but output file was not produced or is 0 bytes: '$abs_output'.\nStderr:\n$err_text\nStdout:\n$out_text")
        end
    finally
        rm(temp_json; force=true)
        rm(temp_script; force=true)
    end
    return abs_output
end

# -----------------------------------------------------------------------------
# 7. Pure Julia 3D Exporters
# -----------------------------------------------------------------------------

"""
    export_obj_pure(obj::BlenderProteinPayload, path::String)

Exports 3D vertex coordinates and octahedron sphere mesh faces as a Wavefront `.obj` file in pure Julia.
"""
function export_obj_pure(obj::BlenderProteinPayload, path::String)
    abs_path = abspath(path)
    open(abs_path, "w") do f
        println(f, "# Wavefront OBJ generated by BioToolkit.jl export_obj_pure")
        n_verts = size(obj.coordinates, 1)
        r = 0.5
        v_idx = 1
        for i in 1:n_verts
            cx, cy, cz = obj.coordinates[i, 1], obj.coordinates[i, 2], obj.coordinates[i, 3]
            println(f, "g Atom_$i")
            @printf(f, "v %.6f %.6f %.6f\n", cx + r, cy, cz)
            @printf(f, "v %.6f %.6f %.6f\n", cx - r, cy, cz)
            @printf(f, "v %.6f %.6f %.6f\n", cx, cy + r, cz)
            @printf(f, "v %.6f %.6f %.6f\n", cx, cy - r, cz)
            @printf(f, "v %.6f %.6f %.6f\n", cx, cy, cz + r)
            @printf(f, "v %.6f %.6f %.6f\n", cx, cy, cz - r)
            
            println(f, "f $(v_idx) $(v_idx+2) $(v_idx+4)")
            println(f, "f $(v_idx) $(v_idx+3) $(v_idx+4)")
            println(f, "f $(v_idx+1) $(v_idx+2) $(v_idx+4)")
            println(f, "f $(v_idx+1) $(v_idx+3) $(v_idx+4)")
            println(f, "f $(v_idx) $(v_idx+2) $(v_idx+5)")
            println(f, "f $(v_idx) $(v_idx+3) $(v_idx+5)")
            println(f, "f $(v_idx+1) $(v_idx+2) $(v_idx+5)")
            println(f, "f $(v_idx+1) $(v_idx+3) $(v_idx+5)")
            v_idx += 6
        end
    end
    return abs_path
end

"""
    export_gltf_pure(obj::AbstractBlenderObject, path::String)

Exports 3D payload scene graph into binary glTF/GLB format using Blender headless rendering if available.
"""
function export_gltf_pure(obj::AbstractBlenderObject, path::String)
    abs_path = abspath(path)
    if is_blender_available()
        scene = BlenderScene(objects=[obj])
        return render_blender_headless(scene, abs_path; format=:glb)
    else
        error("glTF/GLB export requires Blender executable on PATH or ENV[\"BIOTOOLKIT_BLENDER_PATH\"]. Please install Blender or use export_obj_pure.")
    end
end

# -----------------------------------------------------------------------------
# 8. High-Level Helper API
# -----------------------------------------------------------------------------

"""
    extract_coordinates(obj::AbstractBlenderObject)

Extracts N x 3 matrix of spatial coordinates from any Blender payload object.
Throws an error if not implemented for a specific custom payload.
"""
function extract_coordinates end

extract_coordinates(obj::AbstractBlenderObject) = error("extract_coordinates is not implemented for payload type $(typeof(obj))")

extract_coordinates(obj::BlenderProteinPayload) = obj.coordinates
extract_coordinates(obj::BlenderSpatialPayload) = obj.coordinates
extract_coordinates(obj::BlenderContactPayload) = obj.coordinates
extract_coordinates(obj::BlenderHiCPayload) = obj.spline_points
extract_coordinates(obj::BlenderRNAVelocityPayload) = obj.origins
extract_coordinates(obj::BlenderTextPayload) = reshape(Float64[obj.location[1], obj.location[2], obj.location[3]], 1, 3)

function extract_coordinates(obj::BlenderGWASPayload)
    nr, nc = size(obj.height_matrix)
    pts = zeros(Float64, nr * nc, 3)
    idx = 1
    for r in 1:nr, c in 1:nc
        pts[idx, 1] = Float64(r)
        pts[idx, 2] = Float64(c)
        pts[idx, 3] = Float64(obj.height_matrix[r, c])
        idx += 1
    end
    return pts
end

function extract_coordinates(obj::BlenderHMMPayload)
    if size(obj.posterior_matrix, 1) > 0 && size(obj.posterior_matrix, 2) > 0
        nr, nc = size(obj.posterior_matrix)
        pts = zeros(Float64, nr * nc, 3)
        idx = 1
        for r in 1:nr, c in 1:nc
            pts[idx, 1] = Float64(r)
            pts[idx, 2] = Float64(c)
            pts[idx, 3] = Float64(obj.posterior_matrix[r, c] * 10.0)
            idx += 1
        end
        return pts
    else
        n_states = length(obj.states)
        pts = zeros(Float64, max(n_states, 1), 3)
        radius = max(Float64(n_states) * 1.5, 5.0)
        for i in 1:n_states
            angle = 2.0 * pi * (i - 1) / max(n_states, 1)
            pts[i, 1] = radius * cos(angle)
            pts[i, 2] = radius * sin(angle)
            pts[i, 3] = 0.0
        end
        return pts
    end
end

"""
    auto_frame_camera(objects::Vector{<:AbstractBlenderObject}; fov=45.0)

Automatically calculates the 3D bounding box and centroid of scene objects to position and orient the camera.
"""
function auto_frame_camera(objects::Vector{<:AbstractBlenderObject}; fov::Float64=45.0)
    all_pts = Matrix{Float64}(undef, 0, 3)
    for obj in objects
        pts = extract_coordinates(obj)
        all_pts = vcat(all_pts, pts)
    end
    if size(all_pts, 1) == 0
        return BlenderCamera(location=(0.0, -40.0, 10.0), rotation=(60.0, 0.0, 0.0), fov=fov)
    end

    cx = mean(all_pts[:, 1])
    cy = mean(all_pts[:, 2])
    cz = mean(all_pts[:, 3])

    dists = sqrt.((all_pts[:, 1] .- cx).^2 .+ (all_pts[:, 2] .- cy).^2 .+ (all_pts[:, 3] .- cz).^2)
    max_r = max(maximum(dists), 1.0)

    cam_y = cy - 2.5 * max_r
    cam_z = cz + max_r

    dy = cy - cam_y
    dz = cz - cam_z
    rot_x = 90.0 - rad2deg(atan(dz, dy))

    return BlenderCamera(location=(cx, cam_y, cam_z), rotation=(rot_x, 0.0, 0.0), fov=fov)
end

auto_frame_camera(obj::AbstractBlenderObject; fov::Float64=45.0) = auto_frame_camera([obj]; fov=fov)

"""
    quick_render(domain_obj; output_path="render.png", format=:png, engine=:cycles, samples=128)

Renders any BioToolkit domain object or payload in a single line. Automatically frames camera and builds scene.
"""
function quick_render(domain_obj; output_path::String="render.png", format::Symbol=:png, engine::Symbol=:cycles, samples::Int=128)
    payload = domain_obj isa AbstractBlenderObject ? domain_obj : to_blender_payload(domain_obj)
    cam = auto_frame_camera(payload)
    scene = BlenderScene(objects=[payload], camera=cam, engine=engine, samples=samples)
    return render_blender_headless(scene, output_path; format=format)
end

"""
    open_in_blender(domain_obj_or_scene)

Launches Blender in GUI interactive mode pre-populated with the generated 3D scene.
Writes temporary files to a documented directory (`tempdir()/biotoolkit_blender`).
"""
function open_in_blender(domain_obj_or_scene)
    blender_bin = detect_blender_executable()
    if blender_bin === nothing
        error("Blender executable not found. Install Blender or set ENV[\"BIOTOOLKIT_BLENDER_PATH\"].")
    end

    scene = if domain_obj_or_scene isa BlenderScene
        domain_obj_or_scene
    elseif domain_obj_or_scene isa AbstractBlenderObject
        payload = domain_obj_or_scene
        cam = auto_frame_camera(payload)
        BlenderScene(objects=[payload], camera=cam)
    else
        payload = to_blender_payload(domain_obj_or_scene)
        cam = auto_frame_camera(payload)
        BlenderScene(objects=[payload], camera=cam)
    end

    tmp_dir = joinpath(tempdir(), "biotoolkit_blender")
    mkpath(tmp_dir)
    stamp = string(round(Int, time() * 1000))
    temp_json = joinpath(tmp_dir, "scene_$stamp.json")
    write_payload_json(scene, temp_json)
    script_str = build_bpy_script(scene; data_path=temp_json)
    temp_script = joinpath(tmp_dir, "bootstrap_$stamp.py")
    write(temp_script, script_str)

    cmd = `$blender_bin --python $temp_script`
    run(cmd, wait=false)
    return true
end

"""
    quick_export_3d(domain_obj, output_path="model.gltf")

Exports any BioToolkit domain object directly into 3D mesh formats (.gltf / .glb / .obj).
"""
function quick_export_3d(domain_obj, output_path::String="model.gltf")
    payload = domain_obj isa AbstractBlenderObject ? domain_obj : to_blender_payload(domain_obj)
    if endswith(lowercase(output_path), ".obj") && payload isa BlenderProteinPayload
        return export_obj_pure(payload, output_path)
    else
        return export_gltf_pure(payload, output_path)
    end
end

"""
    render_turntable(domain_obj_or_scene; output_path="turntable.mp4", frames=36, fps=24, engine=:cycles, samples=64)

Generates a 360-degree turntable camera animation sequence or MP4 video for any BioToolkit object.
"""
function render_turntable(domain_obj_or_scene; output_path::String="turntable.mp4", frames::Int=36, fps::Int=24, engine::Symbol=:cycles, samples::Int=64)
    scene = if domain_obj_or_scene isa BlenderScene
        domain_obj_or_scene
    elseif domain_obj_or_scene isa AbstractBlenderObject
        payload = domain_obj_or_scene
        cam = auto_frame_camera(payload)
        BlenderScene(objects=[payload], camera=cam, engine=engine, samples=samples)
    else
        payload = to_blender_payload(domain_obj_or_scene)
        cam = auto_frame_camera(payload)
        BlenderScene(objects=[payload], camera=cam, engine=engine, samples=samples)
    end

    fmt = endswith(lowercase(output_path), ".mp4") ? :mp4 : (endswith(lowercase(output_path), ".gif") ? :gif : :png)
    return render_blender_headless(scene, output_path; format=fmt, turntable=true, frames=frames, fps=fps)
end

end # module BlenderIntegrator
