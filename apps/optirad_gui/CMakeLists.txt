add_executable(optirad_gui_app main.cpp)

target_link_libraries(optirad_gui_app PRIVATE
    optirad_gui
)

# Set output name to optirad_gui (without _app suffix)
set_target_properties(optirad_gui_app PROPERTIES OUTPUT_NAME optirad_gui)
