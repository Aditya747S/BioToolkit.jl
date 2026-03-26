module BioToolkitCUDAExt

using BioToolkit
using CUDA

function __init__()
    @eval BioToolkit const CUDA = $CUDA
end

end