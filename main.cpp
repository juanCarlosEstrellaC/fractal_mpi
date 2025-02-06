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

uint32_t _bswap32(uint32_t a) {
    return
    ((a & 0X000000FF) << 24) |  // Mueve el byte menos significativo al más significativo
        ((a & 0X0000FF00) <<  8) |  // Mueve el segundo byte a la posición correcta
            ((a & 0x00FF0000) >>  8) |  // Mueve el tercer byte hacia la derecha
                ((a & 0xFF000000) >> 24);  // Mueve el byte más significativo al menos significativo
}

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

void init() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW!" << std::endl;
        exit(-1);
    }

    window = glfwCreateWindow(WIDTH, HEIGHT, "Fractal MPI Grupal", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        std::cerr << "Failed to create GLFW window!" << std::endl;
        exit(-1);
    }

    glfwSetKeyCallback(window, [](GLFWwindow* window, auto key, int scancode, int action, int mods) {
        if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
            glfwSetWindowShouldClose(window, GLFW_TRUE);  // Cierra la ventana al presionar ESC
    });

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

    //glfwMakeContextCurrent(window);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-1,1,-1,1,-1,1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glfwSwapInterval(1);
    initTextures();
    glEnable(GL_TEXTURE_2D);
}

// Función para calcular Mandelbrot
int divergente(double cx, double cy) {
    int iter = 0;
    double vx = cx, vy = cy;
    while (vx * vx + vy * vy <= 4.0 && iter < max_iterations) {
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

void mandelbrotCpu(int rank, int nprocs, std::vector<int>& buffer_local) {
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
}

void aux(int argc, char* argv[]){
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    int filas_por_rank = HEIGHT / nprocs;
    int altura_local = (rank == nprocs - 1) ? (HEIGHT - filas_por_rank * (nprocs - 1)) : filas_por_rank;
    std::vector<int> buffer_local(altura_local * WIDTH);

    if (rank == 0) {
        pixel_buffer = new GLuint[TOTAL];
    }

    mandelbrotCpu(rank, nprocs, buffer_local);

    std::vector<int> recvcounts(nprocs);
    std::vector<int> displs(nprocs);
    for (int i = 0; i < nprocs; ++i) {
        recvcounts[i] = (i == nprocs - 1) ? (HEIGHT - filas_por_rank * (nprocs - 1)) * WIDTH : filas_por_rank * WIDTH;
        displs[i] = i * filas_por_rank * WIDTH;
    }

    MPI_Gatherv(buffer_local.data(), buffer_local.size(), MPI_INT,
                pixel_buffer, recvcounts.data(), displs.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    /*if (rank == 0) {
        run();
        delete[] pixel_buffer;
    }*/

    MPI_Finalize();
}

void paint(int argc, char* argv[]) {
    fps.update();

    aux(argc,argv);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixel_buffer);

    glBegin(GL_QUADS);{
        glTexCoord2f(0, 1);
        glVertex3f(-1, -1, 0);

        glTexCoord2f(0, 0);
        glVertex3f(-1, 1, 0);

        glTexCoord2f(1, 0);
        glVertex3f(1, 1, 0);

        glTexCoord2f(1, 1);
        glVertex3f(1, -1, 0);
    }
    glEnd();
}

void loop(int argc, char* argv[]) {
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        paint(argc, argv);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
}

void run(int argc, char* argv[]) {
    init();
    loop( argc,  argv);
    glfwTerminate();
}


int main(int argc, char* argv[]) {
    run( argc,  argv);
    delete[] pixel_buffer;
    return 0;
}
