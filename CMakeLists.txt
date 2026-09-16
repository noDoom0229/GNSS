cmake_minimum_required(VERSION 3.15)
project(GNSS_RTK CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SOURCES
    # 原有 SPP 工程
    src/1MatrixAndVector.cpp
    src/2TimeSystemTrans.cpp
    src/3CoordSystemTrans.cpp
    src/21BinaryDataDecode.cpp
    src/32SatellitePVT.cpp
    src/33ErrorCorrect.cpp
    src/34SPP35SPV.cpp
    src/36OutPutResult.cpp
    # RTK 新增
    src/config.cpp
    src/RTK.cpp
    src/Kalman.cpp
    src/lambdaN.cpp
    src/writeObs.cpp
    src/writeToFile.cpp
    src/main.cpp
)

add_executable(GNSS_RTK ${SOURCES})
target_include_directories(GNSS_RTK PRIVATE include)

if(MSVC)
    # 源文件为 UTF-8；屏蔽 fopen/inet_addr 等安全警告
    target_compile_options(GNSS_RTK PRIVATE /utf-8 /W3)
    target_compile_definitions(GNSS_RTK PRIVATE _CRT_SECURE_NO_WARNINGS _WINSOCK_DEPRECATED_NO_WARNINGS)
endif()

if(WIN32)
    target_link_libraries(GNSS_RTK PRIVATE ws2_32)
endif()

# 调试运行时以工程根目录为工作目录，使 config/data/output 相对路径可用
set_property(TARGET GNSS_RTK PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")
