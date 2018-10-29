# Additional targets to perform clang-format/clang-tidy
# Get all project files
file(GLOB_RECURSE ALL_SOURCE_FILES
     *.cpp *.h *.cxx *.hxx *.hpp *.cc
     )

set(CLANG_FORMAT_EXCLUDE_PATTERNS ${CLANG_FORMAT_EXCLUDE_PATTERNS} "/CMakeFiles/" "cmake")

# get all project files file
foreach (SOURCE_FILE ${ALL_SOURCE_FILES})
    foreach (EXCLUDE_PATTERN ${CLANG_FORMAT_EXCLUDE_PATTERNS})
        string(FIND ${SOURCE_FILE} ${EXCLUDE_PATTERN} EXCLUDE_FOUND)
        if (NOT ${EXCLUDE_FOUND} EQUAL -1)
            list(REMOVE_ITEM ALL_SOURCE_FILES ${SOURCE_FILE})
        endif ()
    endforeach ()
endforeach ()

# Adding clang-format target if executable is found
find_program(CLANG_FORMAT "clang-format")
if(CLANG_FORMAT)
  add_custom_target(
    clang-format
    COMMAND ${CLANG_FORMAT} -style=file
    -i
    ${ALL_SOURCE_FILES}
    )

    add_custom_target(
        check-format
    )

    foreach(SOURCE_FILE ${ALL_SOURCE_FILES})
      add_custom_command(TARGET check-format PRE_BUILD
                         COMMAND
                             ${CLANG_FORMAT} -style=file ${SOURCE_FILE} | diff -U5 ${SOURCE_FILE} -
                         )
    endforeach()
endif()