#!/opt/julia0.6/julia

include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFTypes.jl")
include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFFileprc.jl")
include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFSupCalc.jl")
include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFTbCalc.jl")
include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFHamiCalc.jl")
include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFUtils.jl")
# include("src/DFGraphics.jl")
include("/home/ponet/Documents/PhD/Coding/Julia/DFTools.jl/src/DFVisualize.jl")

using GLVisualize, GLWindow, GeometryTypes, Colors, GLAbstraction,GLVisualize.GLVisualizeShader

x = create_TB_model("/home/ponet/Documents/PhD/GeTe/NSOC/paperxsf/","/home/ponet/Documents/PhD/GeTe/SOC/GeTe_bands.out",[[PhysAtom(0.0,0.0,-0.0239129,-0.155854) for i=1:4]...,[PhysAtom(0.0,0.0,5.5540692,0.318205) for i=1:4]...],Float64);

visualize_wfc(x.wfcs[1],0.35)

# using Meshing,GeometryTypes,ColorTypes,GLVisualize

# sdf = SignedDistanceField(HyperRectangle(Vec(-1,-1,-1.),Vec(2,2,2.))) do v
#     sqrt(sum(dot(v,v))) - 1 # sphere
# end;

# # meshing marching_cubes test 
# m = marching_cubes(sdf,0.1);
# it_colors = [[RGBA(1.0,0.0,0.0,0.6) for i=1:3];[RGBA(0.0,0.0,1.0,0.6) for i=1:3]];
# m_test = GLNormalVertexcolorMesh(vertices=m.vertices,faces=m.faces,color=[it_colors[(i-1)%6+1] for i = 1:length(m.vertices)]);

# screen = glscreen()
# _view(visualize(m_test),camera=:perspective)
# renderloop(screen)