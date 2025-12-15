cmake_minimum_required(VERSION 3.10)
project(SupernovaSim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- КРИТИЧЕСКИ ВАЖНО ДЛЯ RYZEN 7 ---
# -O3: Включает агрессивную оптимизацию циклов
# -march=native: Разрешает использовать AVX-512 и другие инструкции твоего CPU
# -ffast-math: Упрощает математику (может чуть снизить точность, но ускорит физику)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -ffast-math")

find_package(SFML 2.5 COMPONENTS graphics window system REQUIRED)

add_executable(main src/main.cpp)
target_link_libraries(main sfml-graphics sfml-window sfml-system)