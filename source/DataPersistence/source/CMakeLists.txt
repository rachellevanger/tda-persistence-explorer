# This CMakeLists.txt file is called recursively from the one in ./..

# add the executable

set ( LIBS      ${LIBS}
                ${X11_LIBRARIES}
                ${CMAKE_THREAD_LIBS_INIT}
                ${Boost_LIBRARIES} )

add_executable ( DataPersistence DataPersistence.cpp )
target_link_libraries ( DataPersistence ${LIBS} )

install(TARGETS DataPersistence
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)