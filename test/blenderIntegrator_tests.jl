using Test
using BioToolkit
using SparseArrays
import Base.Filesystem: tempname

@testset "BlenderIntegrator Core Tests" begin
    @testset "Scene Graph Data Structures & Type Hierarchy" begin
        cam = BlenderCamera(location=(0.0, -40.0, 10.0), rotation=(45.0, 0.0, 0.0), fov=50.0)
        @test cam.location == (0.0, -40.0, 10.0)
        @test cam.rotation == (45.0, 0.0, 0.0)
        @test cam.fov == 50.0

        light = BlenderLight(:sun, (10.0, 10.0, 20.0), 500.0)
        @test light.light_type == :sun
        @test light.energy == 500.0

        mat = BlenderMaterial(name="GlassMat", color=(0.1, 0.5, 0.9, 0.8), roughness=0.1, transmission=0.9)
        @test mat.name == "GlassMat"
        @test mat.transmission == 0.9

        scene = BlenderScene(camera=cam, lights=[light], engine=:cycles, samples=64)
        @test scene.engine == :cycles
        @test scene.samples == 64
        @test length(scene.lights) == 1

        # Type Hierarchy test
        coords = [0.0 0.0 0.0; 1.0 1.0 1.0]
        prot_pl = BlenderProteinPayload("Prot", coords, ["C", "C"], [1, 2], :vdw, :element, mat)
        @test prot_pl isa AbstractBlenderPayload
        @test prot_pl isa AbstractBlenderObject
    end

    @testset "Domain Payloads & Generic Hooks" begin
        atom1 = Atom(1, "CA", 0.0, 0.0, 0.0; element="C")
        atom2 = Atom(2, "CA", 1.5, 2.5, 3.5; element="C")
        res = Residue("ALA", 1, ' ', [atom1, atom2])
        chain = Chain("A", [res])
        model = Model(1, [chain])
        struct_obj = Structure("1TEST", [model], Dict{String,String}())

        payload = to_blender_payload(struct_obj; style=:cartoon)
        @test payload isa BlenderProteinPayload
        @test payload.name == "1TEST"
        @test size(payload.coordinates, 1) == 2

        # Contact Payload
        coords = [0.0 0.0 0.0; 5.0 5.0 5.0; 10.0 10.0 10.0]
        contact_payload = BlenderContactPayload("DCA_Contacts", coords, [(1, 2), (2, 3)], [0.95, 0.88], 0.2, BlenderMaterial())
        @test contact_payload.name == "DCA_Contacts"
        @test length(contact_payload.contact_pairs) == 2

        # Hi-C Payload
        hic_payload = BlenderHiCPayload("Chromatin_3D", coords, [0.1, -0.5, 0.8], [(1, 3)], 0.5, BlenderMaterial())
        @test hic_payload.name == "Chromatin_3D"
        @test size(hic_payload.spline_points, 1) == 3

        # RNA Velocity Payload
        velo_payload = BlenderRNAVelocityPayload("RNA_Velo", coords, coords .* 0.1, 1.0, BlenderMaterial())
        @test velo_payload.name == "RNA_Velo"
        @test size(velo_payload.origins, 1) == 3

        # GWAS Payload
        gwas_h = [0.0 1.5; 2.5 5.0]
        gwas_payload = BlenderGWASPayload("GWAS_Terrain", gwas_h, ["chr1", "chr2"], ["traitA", "traitB"], BlenderMaterial())
        @test gwas_payload.name == "GWAS_Terrain"
        @test size(gwas_payload.height_matrix) == (2, 2)

        # HMM Payload
        hmm_trans = [0.9 0.1; 0.2 0.8]
        hmm_post = [0.8 0.2; 0.3 0.7]
        hmm_payload = BlenderHMMPayload("HMM_Model", ["S1", "S2"], hmm_trans, hmm_post, BlenderMaterial())
        @test hmm_payload.name == "HMM_Model"
        @test size(hmm_payload.transition_matrix) == (2, 2)

        # to_blender_scene method test
        scene_from_obj = to_blender_scene(contact_payload; engine=:eevee)
        @test scene_from_obj isa BlenderScene
        @test scene_from_obj.engine == :eevee

        # Domain to_blender_payload dispatch tests
        cmap = ContactMap([0.0 0.9; 0.9 0.0], [1, 2])
        cmap_pl = to_blender_payload(cmap, coords)
        @test cmap_pl isa BlenderContactPayload

        hic_mat = HiCContactMatrix(SparseArrays.sparse([1.0 0.5; 0.5 1.0]), "chr1", GenomicInterval[], 10000, true, "ICE", Dict{String,Any}())
        hic_pl = to_blender_payload(hic_mat)
        @test hic_pl isa BlenderHiCPayload

        gwas_res = GWASResult(["rs1"], ["chr1"], [100], [("A", "G")], ["GENE1"], [0.5], [0.1], [5.0], [1e-6], 1000, ["cov1"], "trait1", "linear")
        gwas_pl = to_blender_payload(gwas_res)
        @test gwas_pl isa BlenderGWASPayload
    end

    @testset "Material Serialization Consistency & JSON Formatting" begin
        mat = BlenderMaterial(name="SpecMat", color=(0.2, 0.4, 0.6, 1.0), metallic=0.5, roughness=0.3, transmission=0.1, emission_energy=0.2)
        coords = [0.0 0.0 0.0; 1.0 1.0 1.0]
        
        spatial_pl = BlenderSpatialPayload("Spatial", coords, ["C1", "C2"], [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0)], 0.5, mat)
        spatial_dict = to_python_dict(spatial_pl)
        @test haskey(spatial_dict["material"], "metallic")
        @test spatial_dict["material"]["metallic"] == 0.5
        @test spatial_dict["material"]["transmission"] == 0.1

        contact_pl = BlenderContactPayload("Contact", coords, [(1, 2)], [0.9], 0.2, mat)
        contact_dict = to_python_dict(contact_pl)
        @test haskey(contact_dict["material"], "transmission")
    end

    @testset "Extract Coordinates & Camera Math" begin
        coords = [10.0 0.0 0.0; -10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 -10.0 0.0]
        mat = BlenderMaterial(name="TestMat")
        prot_payload = BlenderProteinPayload("ProtCentroid", coords, ["C", "C", "C", "C"], [1, 2, 3, 4], :cartoon, :element, mat)

        @test extract_coordinates(prot_payload) == coords
        
        # Unsupported custom payload should throw error
        struct CustomPayload <: AbstractBlenderObject end
        @test_throws ErrorException extract_coordinates(CustomPayload())

        cam = auto_frame_camera(prot_payload; fov=55.0)
        @test cam isa BlenderCamera
        @test cam.fov == 55.0
        @test cam.location[1] == 0.0
        @test cam.location[2] < 0.0
    end

    @testset "Bootstrap + JSON Payload Serialization" begin
        coords = [1.0 2.0 3.0; 4.0 5.0 6.0]
        mat = BlenderMaterial(name="TestMat")
        prot_payload = BlenderProteinPayload("ProtA", coords, ["C", "N"], [1, 2], :cartoon, :element, mat)
        scene = BlenderScene(objects=[prot_payload])

        tmp_json = tempname() * ".json"
        try
            write_payload_json(scene, tmp_json)
            @test isfile(tmp_json)
            json_content = read(tmp_json, String)
            @test occursin("ProtA", json_content)
            @test occursin("protein", json_content)
            @test occursin("\"fov\":45", json_content)
        finally
            rm(tmp_json; force=true)
        end
    end

    @testset "Pure Julia Script Generation & Features" begin
        coords = [1.0 2.0 3.0; 4.0 5.0 6.0]
        mat = BlenderMaterial(name="TestMat")
        prot_payload = BlenderProteinPayload("ProtA", coords, ["C", "N"], [1, 2], :cartoon, :element, mat)
        gwas_h = [0.0 1.5; 2.5 5.0]
        gwas_payload = BlenderGWASPayload("GWAS_Terrain", gwas_h, ["chr1", "chr2"], ["traitA", "traitB"], mat)
        hic_payload = BlenderHiCPayload("HiC", coords, [0.1, 0.2], [(1, 2)], 0.2, mat)

        scene = BlenderScene(objects=[prot_payload, gwas_payload, hic_payload])
        script_str = build_bpy_script(scene; data_path="payload.json", output_path="test_render.png", format=:png)

        @test occursin("def apply_material", script_str)
        @test occursin("gwas_terrain", script_str)
        @test occursin("spline.use_endpoint_u = True", script_str)
        @test occursin("color_scheme == 'element'", script_str)
        @test occursin("CYCLES", script_str)

        # Output path JSON escaping check
        @test occursin("test_render.png", script_str)

        # Write script to disk test
        tmp_script = tempname() * ".py"
        try
            write_blender_script(scene, tmp_script; data_path="payload.json")
            @test isfile(tmp_script)
            @test filesize(tmp_script) > 0
        finally
            rm(tmp_script; force=true)
        end
    end

    @testset "Protein Color Schemes Script Generation" begin
        coords = [0.0 0.0 0.0; 1.0 1.0 1.0; 2.0 2.0 2.0]
        mat = BlenderMaterial(name="ProtMat")
        
        # Test :rainbow
        payload_rainbow = BlenderProteinPayload("RainbowProt", coords, ["C", "C", "C"], [1, 2, 3], :cartoon, :rainbow, mat)
        scene_rainbow = BlenderScene(objects=[payload_rainbow])
        script_rainbow = build_bpy_script(scene_rainbow; data_path="dummy.json")
        @test occursin("color_scheme == 'rainbow'", script_rainbow)

        # Test :secondary_structure
        payload_ss = BlenderProteinPayload("SSProt", coords, ["C", "C", "C"], [1, 2, 3], :cartoon, :secondary_structure, mat)
        scene_ss = BlenderScene(objects=[payload_ss])
        script_ss = build_bpy_script(scene_ss; data_path="dummy.json")
        @test occursin("color_scheme == 'secondary_structure'", script_ss)

        # Test :sasa
        payload_sasa = BlenderProteinPayload("SASAProt", coords, ["C", "C", "C"], [1, 2, 3], :cartoon, :sasa, mat, [10.0, 50.0, 90.0])
        scene_sasa = BlenderScene(objects=[payload_sasa])
        script_sasa = build_bpy_script(scene_sasa; data_path="dummy.json")
        @test occursin("color_scheme == 'sasa'", script_sasa)
    end

    @testset "Contact Pair Spline Radius & HMM State Networks" begin
        coords = [0.0 0.0 0.0; 10.0 0.0 0.0]
        mat = BlenderMaterial(name="DCA_Mat")
        contact_pl = BlenderContactPayload("DCA_Net", coords, [(1, 2)], [0.95], 0.25, mat)
        scene_contact = BlenderScene(objects=[contact_pl])
        script_contact = build_bpy_script(scene_contact; data_path="dummy.json")
        @test occursin("coupling_scores", script_contact)
        @test occursin("spline.points[0].radius = r_val", script_contact)

        # HMM State Networks
        hmm_trans = [0.9 0.1; 0.2 0.8]
        hmm_post = [0.8 0.2; 0.3 0.7]
        hmm_pl = BlenderHMMPayload("HMM_Model", ["S1", "S2"], hmm_trans, hmm_post, mat)
        scene_hmm = BlenderScene(objects=[hmm_pl])
        script_hmm = build_bpy_script(scene_hmm; data_path="dummy.json")
        @test occursin("trans_matrix", script_hmm)
        @test occursin("primitive_uv_sphere_add", script_hmm)
    end

    @testset "RNA Velocity Vector Field Script Generation" begin
        coords = [0.0 0.0 0.0; 1.0 1.0 1.0]
        mat = BlenderMaterial(name="VeloMat")
        velo_pl = BlenderRNAVelocityPayload("RNAVelo", coords, coords .* 0.5, 2.0, mat)
        scene_velo = BlenderScene(objects=[velo_pl])
        script_velo = build_bpy_script(scene_velo; data_path="dummy.json")
        @test occursin("curve.bevel_depth = 0.05 * arrow_scale", script_velo)
        @test occursin("primitive_cone_add", script_velo)
    end

    @testset "Turntable Animation Script Generation" begin
        mat = BlenderMaterial(name="Mat")
        coords = [0.0 0.0 0.0; 1.0 1.0 1.0]
        prot = BlenderProteinPayload("P1", coords, ["C", "C"], [1, 2], :vdw, :element, mat)
        scene = BlenderScene(objects=[prot])
        
        script_tt = build_bpy_script(scene; output_path="turntable.mp4", turntable=true, frames=48, fps=30, format=:mp4)
        @test occursin("Turntable keyframe animation", script_tt)
        @test occursin("scene.frame_end = 48", script_tt)
        @test occursin("scene.render.fps = 30", script_tt)
        @test occursin("keyframe_insert", script_tt)
        @test occursin("FFMPEG", script_tt)
    end

    @testset "Pure Julia 3D Exporters" begin
        coords = [0.0 1.0 2.0; 3.0 4.0 5.0]
        mat = BlenderMaterial(name="TestMat")
        prot_payload = BlenderProteinPayload("ProtA", coords, ["C", "N"], [1, 2], :vdw, :element, mat)

        tmp_obj = tempname() * ".obj"
        tmp_gltf = tempname() * ".gltf"
        try
            export_obj_pure(prot_payload, tmp_obj)
            @test isfile(tmp_obj)
            obj_str = read(tmp_obj, String)
            @test occursin("Wavefront OBJ", obj_str)
            @test occursin("v 0.500000 1.000000 2.000000", obj_str)
            @test occursin("f 1 3 5", obj_str)

            if !is_blender_available()
                @test_throws ErrorException export_gltf_pure(prot_payload, tmp_gltf)
            end
        finally
            rm(tmp_obj; force=true)
            rm(tmp_gltf; force=true)
        end
    end

    @testset "Blender Binary Autodetection" begin
        @test is_blender_available() isa Bool
    end
end
