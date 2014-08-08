project(ratelesscodes)


add_executable(ratelesscodes main.cpp convcode.cpp convcode.h)

install(TARGETS ratelesscodes RUNTIME DESTINATION bin)
target_link_libraries(ratelesscodes itpp)

add_definitions(-std=c++0x -g)