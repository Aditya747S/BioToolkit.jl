module BioToolkitFluxExt

using BioToolkit
using Flux

function __init__()
    @eval BioToolkit const Flux = $Flux
end

end
