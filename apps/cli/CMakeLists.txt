add_executable(samaya_cli main.cpp)
set_target_properties(samaya_cli PROPERTIES OUTPUT_NAME samaya)
target_link_libraries(samaya_cli PRIVATE samaya::samaya samaya_options)
install(TARGETS samaya_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
