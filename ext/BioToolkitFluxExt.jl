module BioToolkitFluxExt

using BioToolkit
using Flux

function __init__()
    BioToolkit._FLUX_MODULE[] = Flux
end

end
