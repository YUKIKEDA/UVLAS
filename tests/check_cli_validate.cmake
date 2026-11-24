cmake_minimum_required(VERSION 3.25)

function(check expected_code diagnostic)
    execute_process(COMMAND "${UVLAS_EXECUTABLE}" validate ${ARGN}
        WORKING_DIRECTORY "${CASE_DIRECTORY}"
        RESULT_VARIABLE code OUTPUT_VARIABLE output ERROR_VARIABLE error
        ENCODING UTF-8 TIMEOUT 5)
    string(REPLACE "\r\n" "\n" output "${output}")
    if(NOT "${code}" STREQUAL "${expected_code}")
        message(FATAL_ERROR "Exit code ${code}, expected ${expected_code}. stderr: ${error}")
    endif()
    if(expected_code EQUAL 0)
        if(NOT output STREQUAL "Case validation succeeded.\n" OR NOT error STREQUAL "")
            message(FATAL_ERROR "Expected only a success message. stdout: ${output} stderr: ${error}")
        endif()
    else()
        string(FIND "${error}" "${diagnostic}" diagnostic_position)
        if(NOT output STREQUAL "" OR diagnostic_position LESS 0)
            message(FATAL_ERROR "Expected ${diagnostic}. stdout: ${output} stderr: ${error}")
        endif()
    endif()
endfunction()

file(MAKE_DIRECTORY "${CASE_DIRECTORY}")
file(READ "${CASE_EXAMPLES}" examples)
string(JSON structural GET "${examples}" schema_examples structural_static)
set(case_name "検証 case.json")
set(case_path "${CASE_DIRECTORY}/${case_name}")
file(WRITE "${case_path}" "${structural}")
check(0 "" "${case_name}")
check(0 "" "${case_path}")
check(2 "Expected one case path.")
check(2 "manifest acquisition" missing.json)

file(WRITE "${case_path}" "{\"schema\":\"uvlas.case/1\",")
check(2 "JSON document" "${case_name}")
file(WRITE "${case_path}" "{\"schema\":\"uvlas.case/1\",\"schema\":\"uvlas.case/1\"}")
check(2 "JSON document" "${case_name}")

string(JSON broken SET "${structural}" model shells 0 geometry "\"99999\"")
file(WRITE "${case_path}" "${broken}")
check(2 "case schema" "${case_name}")

set(source_path "${CASE_DIRECTORY}/observation.csv")
file(WRITE "${source_path}" "abc")
file(SHA256 "${source_path}" digest)
string(JSON reference GET "${examples}" fragments "ComparisonReference.file")
string(JSON reference SET "${reference}" source source_id "\"900\"")
string(JSON reference SET "${reference}" source path "\"observation.csv\"")
string(JSON reference SET "${reference}" source sha256 "\"${digest}\"")
string(JSON reference SET "${reference}" columns "[\"x\",\"y\",\"z\"]")
string(JSON reference SET "${reference}" uncertainty "[0,0,0]")
string(JSON external SET "${structural}" observations 0 comparison reference "${reference}")
file(WRITE "${case_path}" "${external}")
check(2 "external reference dependency validation is unsupported" "${case_name}")
file(WRITE "${source_path}" "xyz")
check(2 "reference SHA-256 mismatch at /observations/0/comparison/reference/source" "${case_name}")
file(REMOVE "${source_path}")
check(2 "reference acquisition at /observations/0/comparison/reference/source" "${case_name}")

set(mesh [=[$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 3 1 3
2 1 0 3
1
2
3
0 0 0
1 0 0
0 1 0
$EndNodes
$Elements
1 1 1 1
2 1 2 1
1 1 2 3
$EndElements
]=])
file(WRITE "${CASE_DIRECTORY}/skin.msh" "${mesh}")
file(WRITE "${CASE_DIRECTORY}/sidecar.json" "{}")
file(WRITE "${CASE_DIRECTORY}/reference.dat" "surface")
file(SHA256 "${CASE_DIRECTORY}/skin.msh" mesh_digest)
file(SHA256 "${CASE_DIRECTORY}/sidecar.json" sidecar_digest)
file(SHA256 "${CASE_DIRECTORY}/reference.dat" geometry_digest)
string(JSON coupled GET "${examples}" schema_examples coupled_transient)
string(JSON gmsh GET "${examples}" fragments "Geometry.gmsh")
string(JSON gmsh SET "${gmsh}" mesh source_id "\"300\"")
string(JSON gmsh SET "${gmsh}" mesh path "\"skin.msh\"")
string(JSON gmsh SET "${gmsh}" mesh sha256 "\"${mesh_digest}\"")
string(JSON gmsh SET "${gmsh}" sidecar source_id "\"301\"")
string(JSON gmsh SET "${gmsh}" sidecar path "\"sidecar.json\"")
string(JSON gmsh SET "${gmsh}" sidecar sha256 "\"${sidecar_digest}\"")
string(JSON gmsh SET "${gmsh}" reference file source_id "\"302\"")
string(JSON gmsh SET "${gmsh}" reference file path "\"reference.dat\"")
string(JSON gmsh SET "${gmsh}" reference file sha256 "\"${geometry_digest}\"")
string(JSON external SET "${coupled}" model geometry 1 "${gmsh}")
file(WRITE "${case_path}" "${external}")
check(2 "external reference dependency validation is unsupported at /model/geometry/1/sidecar"
    "${case_name}")

string(REPLACE "1 1 2 3" "1 1 2 9" bad_mesh "${mesh}")
file(WRITE "${CASE_DIRECTORY}/skin.msh" "${bad_mesh}")
file(SHA256 "${CASE_DIRECTORY}/skin.msh" mesh_digest)
string(JSON external SET "${external}" model geometry 1 mesh sha256 "\"${mesh_digest}\"")
file(WRITE "${case_path}" "${external}")
file(REMOVE "${CASE_DIRECTORY}/sidecar.json")
check(2 "Gmsh mesh format at byte" "${case_name}")
file(REMOVE "${CASE_DIRECTORY}/skin.msh" "${CASE_DIRECTORY}/reference.dat")

set(deck_path "${CASE_DIRECTORY}/skin.bdf")
set(deck "BEGIN BULK\nGRID,1,,3.0,4.0,5.0\nENDDATA\n")
file(CONFIGURE OUTPUT "${deck_path}" CONTENT "${deck}" @ONLY NEWLINE_STYLE LF)
file(WRITE "${CASE_DIRECTORY}/sidecar.json" "{}")
file(SHA256 "${deck_path}" deck_digest)
string(JSON nastran GET "${examples}" fragments "Geometry.nastran")
string(JSON nastran SET "${nastran}" deck source_id "\"300\"")
string(JSON nastran SET "${nastran}" deck path "\"skin.bdf\"")
string(JSON nastran SET "${nastran}" deck sha256 "\"${deck_digest}\"")
string(JSON nastran SET "${nastran}" sidecar source_id "\"301\"")
string(JSON nastran SET "${nastran}" sidecar path "\"sidecar.json\"")
string(JSON nastran SET "${nastran}" sidecar sha256 "\"${sidecar_digest}\"")
string(JSON external SET "${coupled}" model geometry 1 "${nastran}")
file(WRITE "${case_path}" "${external}")
check(2 "external reference dependency validation is unsupported at /model/geometry/1/sidecar"
    "${case_name}")

string(REPLACE "ENDDATA" "INCLUDE 'child.bdf'\nENDDATA" include_deck "${deck}")
file(CONFIGURE OUTPUT "${deck_path}" CONTENT "${include_deck}" @ONLY NEWLINE_STYLE LF)
check(2 "reference SHA-256 mismatch at /model/geometry/1/deck" "${case_name}")
file(SHA256 "${deck_path}" deck_digest)
string(JSON external SET "${external}" model geometry 1 deck sha256 "\"${deck_digest}\"")
file(WRITE "${case_path}" "${external}")
file(REMOVE "${CASE_DIRECTORY}/sidecar.json")
check(2 "Nastran Bulk Data format at byte 31 at /model/geometry/1/deck" "${case_name}")
file(REMOVE "${deck_path}")

math(EXPR oversized "8 * 1024 * 1024 + 1")
string(REPEAT " " ${oversized} invalid)
file(WRITE "${case_path}" "${invalid}")
check(5 "resource limit exceeded at manifest acquisition" "${case_name}")
file(REMOVE "${case_path}")
