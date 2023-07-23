set(common_files
	CahnHilliard.cpp
	Object.cpp
	models/Landau.cpp
	models/SimpleWertheim.cpp
	utils/Delta.cpp
	)
	
add_library(ch_common ${common_files})
	
add_executable(ch_1D main.cpp ${source_files})
target_link_libraries(ch_1D ch_common spdlog::spdlog)
target_compile_definitions(ch_1D PRIVATE DIM=1)

add_executable(ch_2D main.cpp ${source_files})
target_link_libraries(ch_2D ch_common spdlog::spdlog)
target_compile_definitions(ch_2D PRIVATE DIM=2)

add_compile_options(-Wall -Wextra)
