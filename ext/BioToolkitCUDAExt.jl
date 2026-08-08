module BioToolkitCUDAExt

using BioToolkit
using CUDA

function __init__()
    BioToolkit.eval(:(using CUDA))
end

end