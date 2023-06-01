set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

set(toolchainpath C:/build/23.0/x64/sysroots)

set(CMAKE_C_COMPILER ${toolchainpath}/x86_64-w64-mingw32/usr/bin/x86_64-nilrt-linux/x86_64-nilrt-linux-gcc.exe)
set(CMAKE_CXX_COMPILER ${toolchainpath}/x86_64-w64-mingw32/usr/bin/x86_64-nilrt-linux/x86_64-nilrt-linux-g++.exe)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_SYSROOT ${toolchainpath}/core2-64-nilrt-linux)
set(CMAKE_<LANG>_STANDARD_INCLUDE_DIRECTORIES ${toolchainpath}/core2-64-nilrtlinux/usr/include/c++/10.3.0 ${toolchainpath}/core2-64-nilrt-linux/usr/include/c++/10.3.0/x86_64-nilrtlinux)
set(CMAKE_<LANG>_FLAGS "-Wall -fmessage-length=0")
set(CMAKE_<LANG>_FLAGS_DEBUG "-O0 -g3")
set(CMAKE_<LANG>_FLAGS_RELEASE "-O3")

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE BOTH)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE BOTH)





# project specific information
cmake_minimum_required(VERSION 3.19)
project(ActiveNoiseCancellation)





set(EXECUTABLE_OUTPUT_PATH ../run)
set(CMAKE_BUILD_TYPE Release)
set(HEADER_DIRS include lib/Eigen lib lib/boost src lib/yaml-cpp)
set(LIB_DIRS include lib)
set(LIB libnidaqmx.so m libyaml-cpp.a pthread)
set(DAQMXLIBPATH "C:/Program\ Files\ (x86)/National\ Instruments/Shared/ExternalCompilerSupport/C/lib64/gcc")




#find_library(DAQMXLIBPATH NAMES nidaqmx HNTS "C:/Program\ Files\ (x86)/National\ Instruments/Shared/ExternalCompilerSupport/C/lib64/gcc")

# add_executable(ContGen-ExtClk hardware_test/ContGen-ExtClk.cpp)
# target_include_directories(ContGen-ExtClk PUBLIC include)
# target_link_directories(ContGen-ExtClk PUBLIC lib)
# target_link_libraries(ContGen-ExtClk PUBLIC libnidaqmx.so m)

# add_executable(ContGen-ExtClk-Update hardware_test/ContGen-ExtClk-Update.cpp)
# target_include_directories(ContGen-ExtClk-Update PUBLIC include)
# target_link_directories(ContGen-ExtClk-Update PUBLIC lib)
# target_link_libraries(ContGen-ExtClk-Update PUBLIC libnidaqmx.so m)

# add_executable(ConRW hardware_test/ConRW.cpp )
# target_include_directories(ConRW PUBLIC include include lib/Eigen lib lib/boost)
# target_link_directories(ConRW PUBLIC lib )
# target_link_libraries(ConRW PUBLIC libnidaqmx.so m)

# add_executable(SOC hardware_test/second_order_oscillator.cpp)
# target_include_directories(SOC PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(SOC PUBLIC lib)
# target_link_libraries(SOC PUBLIC libnidaqmx.so m)

# add_executable(cppVersion hardware_test/cppVersion.cpp)

# add_executable(taskLoop hardware_test/taskLoop.cpp)
# target_include_directories(taskLoop PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(taskLoop PUBLIC lib)
# target_link_libraries(taskLoop PUBLIC libnidaqmx.so m)

# add_executable(distRej hardware_test/distRej.cpp)
# target_include_directories(distRej PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(distRej PUBLIC lib)
# target_link_libraries(distRej PUBLIC libnidaqmx.so m)

# add_executable(scanPhase hardware_test/scanPhase.cpp)
# target_include_directories(scanPhase PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(scanPhase PUBLIC lib)
# target_link_libraries(scanPhase PUBLIC libnidaqmx.so m)

# add_executable(scanAmp hardware_test/scanAmp.cpp)
# target_include_directories(scanAmp PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(scanAmp PUBLIC lib)
# target_link_libraries(scanAmp PUBLIC libnidaqmx.so m)

# add_executable(distRejTemplate hardware_test/distRejTemplate.cpp)
# target_include_directories(distRejTemplate PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(distRejTemplate PUBLIC lib)
# target_link_libraries(distRejTemplate PUBLIC libnidaqmx.so m)

# add_executable(distRejSwitching hardware_test/distRejSwitching.cpp)
# target_include_directories(distRejSwitching PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(distRejSwitching PUBLIC lib)
# target_link_libraries(distRejSwitching PUBLIC libnidaqmx.so m)

# add_executable(distRejAFC hardware_test/distRejAFC.cpp)
# target_include_directories(distRejAFC PUBLIC include lib/Eigen lib lib/boost)
# target_link_directories(distRejAFC PUBLIC lib)
# target_link_libraries(distRejAFC PUBLIC libnidaqmx.so m)


add_executable(daq_wrapper src/daq_wrapper.cpp src/xjk/controlFunction.cpp)
target_include_directories(daq_wrapper PUBLIC ${HEADER_DIRS} src/xjk)
target_link_directories(daq_wrapper PUBLIC ${LIB_DIRS})
target_link_libraries(daq_wrapper PUBLIC ${LIB})

add_executable(mtcontroller src/daq_wrapper.cpp src/hgq/controlFunction.cpp)
target_include_directories(mtcontroller PUBLIC ${HEADER_DIRS} src/hgq)
target_link_directories(mtcontroller PUBLIC ${LIB_DIRS})
target_link_libraries(mtcontroller PUBLIC ${LIB})

add_executable(mtcontroller_sim src/daq_sim.cpp src/hgq/controlFunction.cpp)
target_include_directories(mtcontroller_sim PUBLIC ${HEADER_DIRS} src/hgq)
target_link_directories(mtcontroller_sim PUBLIC ${LIB_DIRS})
target_link_libraries(mtcontroller_sim PUBLIC ${LIB})

add_executable(daq_sine src/daq_sine.cpp)
target_include_directories(daq_sine PUBLIC ${HEADER_DIRS})
target_link_directories(daq_sine PUBLIC ${LIB_DIRS})
target_link_libraries(daq_sine PUBLIC ${LIB})