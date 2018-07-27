cmake_minimum_required(VERSION 3.12)
project(Voxelyze)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(Utils)

add_executable(Voxelyze
        Utils/Array3D.cpp
        Utils/Array3D.h
        Utils/GL_Utils.cpp
        Utils/GL_Utils.h
        Utils/MarchCube.cpp
        Utils/MarchCube.h
        Utils/Mesh.cpp
        Utils/Mesh.h
        Utils/tinyxml.cpp
        Utils/tinyxml.h
        Utils/tinyxmlerror.cpp
        Utils/tinyxmlparser.cpp
        Utils/Vec3D.h
        Utils/XML_Rip.cpp
        Utils/XML_Rip.h
        Utils/zconf.h
        Utils/zlib.h
        OpenGLInclude.h
        VX_Benchmark.cpp
        VX_Benchmark.h
        VX_Bond.cpp
        VX_Bond.h
        VX_Enums.h
        VX_Environment.cpp
        VX_Environment.h
        VX_FEA.cpp
        VX_FEA.h
        VX_FRegion.cpp
        VX_FRegion.h
        "VX_MeshUtil (copy).cpp"
        VX_MeshUtil.cpp
        VX_MeshUtil.h
        VX_Object.cpp
        VX_Object.h
        VX_Sim.cpp
        VX_Sim.h
        VX_SimGA.cpp
        VX_SimGA.h
        VX_Source.cpp
        VX_Source.h
        VX_Voxel.cpp
        VX_Voxel.h
        VXS_Bond.cpp
        VXS_Bond.h
        VXS_BondCollision.cpp
        VXS_BondCollision.h
        VXS_BondInternal.cpp
        VXS_BondInternal.h
        VXS_SimGLView.cpp
        VXS_SimGLView.h
        VXS_Voxel.cpp
        VXS_Voxel.h)
