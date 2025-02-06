#include "fps_counter.h"
#include <fmt/core.h>
#include <chrono>

fps_counter::fps_counter() {
    frames = 0;
    fps = 0;
    last_time = std::chrono::steady_clock::now(); // Cambiado a steady_clock
}

int fps_counter::get_fps() {
    return this->fps;
}

void fps_counter::update() {
    frames++;

    auto current_time = std::chrono::steady_clock::now(); // Cambiado a steady_clock

    std::chrono::duration<double, std::milli> tiempo = current_time - last_time;

    if (tiempo.count() > 1000) {
        fps = frames;
        frames = 0;
        last_time = current_time;

        fmt::println("FPS: {}", fps);
    }
}
