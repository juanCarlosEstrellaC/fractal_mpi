#include <iostream>
#include <fmt/core.h>
#include <GLFW/glfw3.h>

#include "fps_counter.h"
#include <memory>
#include <iostream>
#include <cmath>
#include <mpi.h>
#include <vector>

uint32_t _bswap32(uint32_t a) {
    return
    ((a & 0X000000FF) << 24) |  // Mueve el byte menos significativo al más significativo
        ((a & 0X0000FF00) <<  8) |  // Mueve el segundo byte a la posición correcta
            ((a & 0x00FF0000) >>  8) |  // Mueve el tercer byte hacia la derecha
                ((a & 0xFF000000) >> 24);  // Mueve el byte más significativo al menos significativo
}

#define WIDTH 800  // Ancho de la ventana
#define HEIGHT 600  // Altura de la ventana
#define TOTAL 480000

// Límites del plano cartesiano para calcular Mandelbrot
const double x_max = 1;
const double x_min = -2;
const double y_min = -1;
const double y_max = 1;

const int max_iterations = 100;  // Iteraciones máximas para determinar convergencia
const int PALETE_SIZE = 16;  // Tamaño de la paleta de colores

const GLuint color_ramp[PALETE_SIZE] = {
    _bswap32(0xFFFF1010),
    _bswap32(0xFFF31017),
    _bswap32(0xFFE8101E),
    _bswap32(0xFFDC1126),
    _bswap32(0xFFD1112D),
    _bswap32(0xFFC51235),
    _bswap32(0xFFBA123C),
    _bswap32(0xFFAE1343),
    _bswap32(0xFFA3134B),
    _bswap32(0xFF971452),
    _bswap32(0xFF8C145A),
    _bswap32(0xFF801461),
    _bswap32(0xFF751568),
    _bswap32(0xFF691570),
    _bswap32(0xFF5E1677),
    _bswap32(0xFF54167D),
};

static GLFWwindow* window = NULL;  // Puntero a la ventana
GLuint textureID;  // ID de textura de OpenGL
GLuint* pixel_buffer = nullptr;  // Buffer para almacenar colores de los píxeles
fps_counter fps;  // Instancia del contador de FPS

void initTextures() {
    glGenTextures(1, &textureID);  // Genera un identificador para la textura
    glBindTexture(GL_TEXTURE_2D, textureID);  // Enlaza la textura como 2D

    // Define una textura vacía con formato RGBA8 y dimensiones WIDTH x HEIGHT
    glTexImage2D(GL_TEXTURE_2D,
            0,
            GL_RGBA8,
            WIDTH, HEIGHT, 0,
            GL_RGBA,
            GL_UNSIGNED_BYTE,
            NULL
    );

    // Configura los filtros de la textura
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glBindTexture(GL_TEXTURE_2D, 0);  // Desenlaza la textura
}

void init() {
    if (!glfwInit()) { // Inicializa GLFW
        std::cerr << "Failed to initialize GLFW!" << std::endl;
        exit(-1);
    }

    // Crea una ventana con las dimensiones especificadas
    window = glfwCreateWindow(WIDTH, HEIGHT, "OpenGl C++", NULL, NULL);
    if (!window) {
        glfwTerminate();
        std::cerr << "Failed to create GLFW window!" << std::endl;
        exit(-1);
    }


    // Configura un callback para manejar eventos de teclado
    glfwSetKeyCallback(window, [](GLFWwindow* window, auto key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GLFW_TRUE);  // Cierra la ventana al presionar ESC
    });

    // Callback para manejar cambios de tamaño de la ventana
    glfwSetFramebufferSizeCallback(window, [](GLFWwindow* window, int width, int height) {
        glViewport(0, 0, width, height);  // Ajusta la vista al nuevo tamaño
    });


    glfwMakeContextCurrent(window);  // Asocia OpenGL a la ventana actual
    std::string version = (char *)glGetString(GL_VERSION);  // Obtiene la versión de OpenGL
    std::string vendor = (char *)glGetString(GL_VENDOR);  // Obtiene el proveedor de GPU
    std::string render = (char *)glGetString(GL_RENDERER);  // Obtiene el nombre del renderizador

    // Imprime información de OpenGL
    fmt::print("OpenGL version supported {}\n", version);
    fmt::print("Vendor : {}\n", vendor);
    fmt::print("Renderer : {}\n", render);

    // Configuración de la proyección ortogonal
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1,1,-1,1,-1,1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_TEXTURE_2D); // Habilita texturas 2D

    glfwSwapInterval(1); // Habilita v-sync (sincronización vertical)

    initTextures(); // Inicializa las texturas
}

int divergente(double cx, double cy) {
    int iter = 0;
    double vx = cx;
    double vy = cy;

    while (iter<max_iterations && (vx * vx + vy * vy) <= 4.0) {
        double tx = vx * vx - vy * vy + cx;
        double ty = 2.0 * vx * vy + cy;
        vx = tx;
        vy = ty;
        iter++;
    }

    if (iter > 0 && iter < max_iterations) {
        int color_idx = iter % PALETE_SIZE;
        return color_ramp[color_idx];
    }

    if ((vx * vx + vy * vy) > 4.0) {
        return color_ramp[0]; // Color para puntos que se escapan
    }

    return 0;  // Converge (fuera del conjunto de Mandelbrot)
}

void mandelbrotCpu() {
    double dx = (x_max - x_min) / WIDTH;
    double dy = (y_max - y_min) / HEIGHT;

    for (int k = 0; k < WIDTH * HEIGHT; k++) {
        int i = k % WIDTH;       // Columna (coordenada x)
        int j = k / WIDTH;       // Fila (coordenada y)

        double x = x_min + i * dx;
        double y = y_max - j * dy;

        pixel_buffer[k] = divergente(x, y);
    }

}

void paint() {


    // Actualizamos el contador de fotogramas por segundo.
    fps.update();

    // Calculamos el conjunto de Mandelbrot y almacenamos los datos en el buffer de píxeles.
    //mandelbrotCpu();

    glEnable(GL_TEXTURE_2D);

    // Enlazamos la textura creada previamente.
    glBindTexture(GL_TEXTURE_2D, textureID);

    // Actualizamos la textura con los colores del buffer de píxeles.
    glTexImage2D(GL_TEXTURE_2D,
            0,                // Nivel base de la textura.
            GL_RGBA,          // Formato interno (RGBA).
            WIDTH, HEIGHT,    // Dimensiones de la textura.
            0,                // Borde (debe ser 0).
            GL_RGBA,          // Formato de los datos.
            GL_UNSIGNED_BYTE, // Tipo de datos en el buffer.
            pixel_buffer);    // Datos del buffer.

    // Dibujamos un cuadrado cubriendo toda la ventana para mostrar la textura.
    glBegin(GL_QUADS);
    {
        // Definimos las coordenadas de textura y los vértices del cuadrado.
        glTexCoord2f(0, 1);
        glVertex3f(-1, -1, 0);

        glTexCoord2f(0, 0);
        glVertex3f(-1, 1, 0);

        glTexCoord2f(1, 0);
        glVertex3f(1, 1, 0);

        glTexCoord2f(1, 1);
        glVertex3f(1, -1, 0);
    }
    glEnd();  // Finalizamos el dibujo del cuadrado.

}

void loop() {
    // Establecemos el color de fondo de la ventana en negro.
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    // Iniciamos el bucle principal de la aplicación.
    while (!glfwWindowShouldClose(window)) {
        // Limpiamos el buffer de color y el buffer de profundidad.
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Llamamos a la función `paint` para calcular y renderizar el conjunto de Mandelbrot.
        paint();

        // Intercambiamos los buffers de la ventana (doble buffering).
        glfwSwapBuffers(window);

        // Procesamos eventos pendientes de la ventana (como entradas del teclado o mouse).
        glfwPollEvents();
    }
}

void run() {
    init();
    loop();
    glfwTerminate();
}


int main(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    std::vector<int> elementos_por_rank(nprocs);     //elementos_por_rank[4] = {120k, 120k, 120k, 120k};
    if (rank == 0){
        pixel_buffer = new GLuint[WIDTH * HEIGHT];

        int e = TOTAL / nprocs;
        int f = TOTAL - e * (nprocs - 1);

        for (int i = 0; i < nprocs; i++){
            elementos_por_rank[i] = (i == nprocs - 1) ? f : e;
            printf("elementos_por_rank[%d] = %d\n", i, elementos_por_rank[i]);
        }
        run();
        delete[] pixel_buffer;
    }

    std::vector<int> buffer_local(elementos_por_rank[rank]);

    std::vector<int> desplazamientos(nprocs);
    desplazamientos[0] = 0;
    for (int i = 1; i < nprocs; i++) {
        desplazamientos[i] = desplazamientos[i - 1] + elementos_por_rank[i - 1];
    }



    double dx = (x_max - x_min) / WIDTH;
    double dy = (y_max - y_min) / HEIGHT;

    int rows_per_process = HEIGHT / nprocs;
    int start_row = rank * rows_per_process;
    int end_row = (rank == nprocs - 1) ? HEIGHT : start_row + rows_per_process;

    for (int j = start_row; j < end_row; ++j) {
        for (int i = 0; i < WIDTH; ++i) {
            double cx = x_min + i * dx;
            double cy = y_max - j * dy;
            buffer_local[(j - start_row) * WIDTH + i] = divergente(cx, cy);
        }
    }



    MPI_Gatherv(buffer_local.data(), elementos_por_rank[rank], MPI_INT,
                pixel_buffer, elementos_por_rank.data(), desplazamientos.data(),
                MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
