# the moldb executable
add_executable(moldb_exe moldb.cpp)
set_target_properties(moldb_exe PROPERTIES OUTPUT_NAME moldb)
target_link_libraries(moldb_exe moldb pthread)
install(TARGETS moldb_exe
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib)

# openbabel service
add_executable(moldb_openbabelservice openbabelservice.cpp)
target_link_libraries(moldb_openbabelservice moldb)
install(TARGETS moldb_openbabelservice
                RUNTIME DESTINATION bin
                LIBRARY DESTINATION lib
                ARCHIVE DESTINATION lib)
