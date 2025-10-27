cmake_minimum_required(VERSION 3.10)
project(D2D)

# 1) 確保 bin 資料夾存在
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# 2) 指定可執行檔輸出到 bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_CXX_STANDARD 17) 

# 確保 result 資料夾存在
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/result)

# 只清空 result 資料夾內容（保留執行檔）
add_custom_target(cleanResult
  COMMENT "清空 result 目錄下的所有檔案，但保留資料夾與執行檔"
  COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_SOURCE_DIR}/result
  COMMAND ${CMAKE_COMMAND} -E make_directory  ${CMAKE_SOURCE_DIR}/result
)

add_custom_target(show
  COMMAND ./bin/ShowResult result/ case/design.rule 
)


# for cgal
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build" FORCE)
endif()

# 尋找 Qt5 Widgets 模組
find_package(Qt5 REQUIRED COMPONENTS Widgets)

# 查找 CGAL 套件
find_package(CGAL REQUIRED COMPONENTS Qt5)

# 自動抓取 main.cpp 以及 src 資料夾底下所有 .cpp
file(GLOB_RECURSE MAIN_SOURCES 
    "${CMAKE_SOURCE_DIR}/src/main.cpp" #1
    # "${CMAKE_SOURCE_DIR}/src/ShowGraph.cpp" #2
    # "${CMAKE_SOURCE_DIR}/src/GenerateCase.cpp" #3
    "${CMAKE_SOURCE_DIR}/libs/*.cpp" 
)

#ShowGraph
#ShowResult  

# file(GLOB_RECURSE GEN_SOURCES
#     "${CMAKE_SOURCE_DIR}/src/BumpGenerator.cpp"
#     "${CMAKE_SOURCE_DIR}/include/*.cpp"
# )
# 添加頭文件目錄
include_directories(${CMAKE_SOURCE_DIR}/include)

# 添加可執行檔案
add_executable(D2D ${MAIN_SOURCES}) #1
# add_executable(ShowResult  ${MAIN_SOURCES}) #2
# add_executable(GenerateCase ${MAIN_SOURCES}) #3

# 連結 Qt5 Widgets 和 CGAL
target_link_libraries(D2D CGAL::CGAL Qt5::Widgets) #1
# target_link_libraries(ShowResult CGAL::CGAL Qt5::Widgets) #2 
# target_link_libraries(GenerateCase CGAL::CGAL Qt5::Widgets) #3