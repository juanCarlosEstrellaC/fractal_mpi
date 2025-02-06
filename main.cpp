#include <iostream>
#include <fmt/core.h>
#include <GLFW/glfw3.h>
#include <vector>
#include <mpi.h>
#include "fps_counter.h"

// Constantes
#define WIDTH 800
#define HEIGHT 600
#define TOTAL (WIDTH * HEIGHT)
const double x_max = 1;
const double x_min = -2;
const double y_min = -1;
const double y_max = 1;
const int max_iterations = 100;
const int PALETE_SIZE = 16;

// Paleta de colores
const GLuint color_ramp[PALETE_SIZE] = {
    0xFFFF1010, 0xFFF31017, 0xFFE8101E, 0xFFDC1126,
    0xFFD1112D, 0xFFC51235, 0xFFBA123C, 0xFFAE1343,
    0xFFA3134B, 0xFF971452, 0xFF8C145A, 0xFF801461,
    0xFF751568, 0xFF691570, 0xFF5E1677, 0xFF54167D
};

// Variables globales
static GLFWwindow* window = nullptr;
GLuint textureID;
GLuint* pixel_buffer = nullptr;
fps_counter fps;


// Inicialización de GLFW y texturas
void initTextures() {
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void initGLFW() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW!" << std::endl;
        exit(EXIT_FAILURE);
    }

    window = glfwCreateWindow(WIDTH, HEIGHT, "Mandelbrot with MPI and GLFW", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        std::cerr << "Failed to create GLFW window!" << std::endl;
        exit(EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    initTextures();
    glEnable(GL_TEXTURE_2D);
}

// Función para calcular Mandelbrot
int calculatePixel(double cx, double cy) {
    double zx = cx, zy = cy;
    int iter = 0;
    while (zx * zx + zy * zy <= 4.0 && iter < max_iterations) {
        double temp = zx * zx - zy * zy + cx;
        zy = 2.0 * zx * zy + cy;
        zx = temp;
        iter++;
    }

    return iter == max_iterations ? 0 : color_ramp[iter % PALETE_SIZE];
}

// Función para renderizar Mandelbrot
void renderMandelbrot(int rank, int nprocs, std::vector<int>& buffer_local) {
    double dx = (x_max - x_min) / WIDTH;
    double dy = (y_max - y_min) / HEIGHT;

    int rows_per_process = HEIGHT / nprocs;
    int start_row = rank * rows_per_process;
    int end_row = (rank == nprocs - 1) ? HEIGHT : start_row + rows_per_process;

    for (int j = start_row; j < end_row; ++j) {
        for (int i = 0; i < WIDTH; ++i) {
            double cx = x_min + i * dx;
            double cy = y_max - j * dy;
            buffer_local[(j - start_row) * WIDTH + i] = calculatePixel(cx, cy);
        }
    }
}

// Función para pintar en pantalla
void paint() {
    fps.update();

    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixel_buffer);

    glBegin(GL_QUADS);
    glTexCoord2f(0, 1); glVertex2f(-1, -1);
    glTexCoord2f(0, 0); glVertex2f(-1, 1);
    glTexCoord2f(1, 0); glVertex2f(1, 1);
    glTexCoord2f(1, 1); glVertex2f(1, -1);
    glEnd();
}

// Bucle principal
void loop() {
    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);
        paint();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

// Programa principal
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int rows_per_process = HEIGHT / nprocs;
    int local_height = (rank == nprocs - 1) ? (HEIGHT - rows_per_process * (nprocs - 1)) : rows_per_process;
    std::vector<int> buffer_local(local_height * WIDTH);

    if (rank == 0) {
        pixel_buffer = new GLuint[TOTAL];
    }

    renderMandelbrot(rank, nprocs, buffer_local);

    std::vector<int> recvcounts(nprocs);
    std::vector<int> displs(nprocs);
    for (int i = 0; i < nprocs; ++i) {
        recvcounts[i] = (i == nprocs - 1) ? (HEIGHT - rows_per_process * (nprocs - 1)) * WIDTH : rows_per_process * WIDTH;
        displs[i] = i * rows_per_process * WIDTH;
    }

    MPI_Gatherv(buffer_local.data(), buffer_local.size(), MPI_INT,
                pixel_buffer, recvcounts.data(), displs.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        initGLFW();
        loop();
        glfwTerminate();
        delete[] pixel_buffer;
    }

    MPI_Finalize();
    return 0;
}
