#ifndef FPS_COUNTER_H
#define FPS_COUNTER_H

#include <chrono> // Incluye std::chrono para la gestión del tiempo.

namespace ch = std::chrono; // Alias para simplificar la escritura.

class fps_counter {
private:
    int frames; // Contador de cuadros procesados.
    int fps;    // Cuadros por segundo calculados.
    ch::time_point<ch::steady_clock> last_time; // Usa steady_clock para tiempos consistentes.

public:
    fps_counter();      // Constructor.
    int get_fps();      // Devuelve el FPS calculado.
    void update();      // Actualiza el contador y calcula el FPS.
};

#endif // FPS_COUNTER_H
