# OptiRad source modules
add_subdirectory(utils)
add_subdirectory(core)
add_subdirectory(geometry)
add_subdirectory(dose)
add_subdirectory(io)
add_subdirectory(optimization)
add_subdirectory(steering)
add_subdirectory(phsp)

if(OPTIRAD_BUILD_GUI)
    add_subdirectory(gui)
endif()
