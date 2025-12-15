#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

// --- КОНСТАНТЫ СИМУЛЯЦИИ ---
const int N = 100; // Размер сетки (теперь 100x100)
const int SIZE = (N + 2) * (N + 2);
const int SCALE = 6; // Размер одного "пикселя" симуляции на экране (100 * 6 = 600px окно)

// Макрос для индексации
int IX(int x, int y) {
    // Ограничиваем координаты, чтобы не вылететь за пределы массива
    if (x < 0) x = 0; if (x > N + 1) x = N + 1;
    if (y < 0) y = 0; if (y > N + 1) y = N + 1;
    return x + (N + 2) * y;
}

// --- КЛАСС ФИЗИКИ (ТОТ ЖЕ, ЧТО БЫЛ, НО С КОНСТАНТАМИ ВНУТРИ) ---
class GasSolver {
public:
    vector<double> u, v, u_prev, v_prev, dens, dens_prev;
    double dt, diff, visc;

    GasSolver() {
        u.resize(SIZE, 0.0); v.resize(SIZE, 0.0);
        u_prev.resize(SIZE, 0.0); v_prev.resize(SIZE, 0.0);
        dens.resize(SIZE, 0.0); dens_prev.resize(SIZE, 0.0);

        dt = 0.1;
        diff = 0.0000; // Диффузия 0, чтобы дым дольше держался
        visc = 0.0000;
    }

    void addDensity(int x, int y, double amount) {
        dens[IX(x, y)] += amount;
    }

    void addVelocity(int x, int y, double amountX, double amountY) {
        u[IX(x, y)] += amountX;
        v[IX(x, y)] += amountY;
    }

    void step() {
        diffuse(1, u_prev, u, visc);
        diffuse(2, v_prev, v, visc);
        project(u_prev, v_prev, u, v);
        advect(1, u, u_prev, u_prev, v_prev);
        advect(2, v, v_prev, u_prev, v_prev);
        project(u, v, u_prev, v_prev);
        diffuse(0, dens_prev, dens, diff);
        advect(0, dens, dens_prev, u, v);
    }

private:
    void set_bnd(int b, vector<double>& x) {
        for (int i = 1; i <= N; i++) {
            x[IX(0, i)] = (b == 1) ? -x[IX(1, i)] : x[IX(1, i)];
            x[IX(N + 1, i)] = (b == 1) ? -x[IX(N, i)] : x[IX(N, i)];
            x[IX(i, 0)] = (b == 2) ? -x[IX(i, 1)] : x[IX(i, 1)];
            x[IX(i, N + 1)] = (b == 2) ? -x[IX(i, N)] : x[IX(i, N)];
        }
        x[IX(0, 0)] = 0.5 * (x[IX(1, 0)] + x[IX(0, 1)]);
        x[IX(0, N + 1)] = 0.5 * (x[IX(1, N + 1)] + x[IX(0, N)]);
        x[IX(N + 1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N + 1, 1)]);
        x[IX(N + 1, N + 1)] = 0.5 * (x[IX(N, N + 1)] + x[IX(N + 1, N)]);
    }

    void diffuse(int b, vector<double>& x, vector<double>& x0, double diff) {
        double a = dt * diff * N * N;
        for (int k = 0; k < 5; k++) { // Уменьшил итерации до 5 для скорости в Real-Time
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] +
                                                       x[IX(i, j - 1)] + x[IX(i, j + 1)])) / (1 + 4 * a);
                }
            }
            set_bnd(b, x);
        }
    }

    void advect(int b, vector<double>& d, vector<double>& d0, vector<double>& u, vector<double>& v) {
        double dt0 = dt * N;
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                double x = i - dt0 * u[IX(i, j)];
                double y = j - dt0 * v[IX(i, j)];
                if (x < 0.5) x = 0.5; if (x > N + 0.5) x = N + 0.5;
                if (y < 0.5) y = 0.5; if (y > N + 0.5) y = N + 0.5;
                int i0 = (int)x; int i1 = i0 + 1;
                int j0 = (int)y; int j1 = j0 + 1;
                double s1 = x - i0; double s0 = 1 - s1;
                double t1 = y - j0; double t0 = 1 - t1;
                d[IX(i, j)] = s0 * (t0 * d0[IX(i0, j0)] + t1 * d0[IX(i0, j1)]) +
                              s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
            }
        }
        set_bnd(b, d);
    }

    void project(vector<double>& u, vector<double>& v, vector<double>& p, vector<double>& div) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                div[IX(i, j)] = -0.5 * (u[IX(i + 1, j)] - u[IX(i - 1, j)] +
                                        v[IX(i, j + 1)] - v[IX(i, j - 1)]) / N;
                p[IX(i, j)] = 0;
            }
        }
        set_bnd(0, div); set_bnd(0, p);
        for (int k = 0; k < 5; k++) { // Уменьшил итерации до 5
            for (int i = 1; i <= N; i++) {
                for (int j = 1; j <= N; j++) {
                    p[IX(i, j)] = (div[IX(i, j)] + p[IX(i - 1, j)] + p[IX(i + 1, j)] +
                                   p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4;
                }
            }
            set_bnd(0, p);
        }
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                u[IX(i, j)] -= 0.5 * N * (p[IX(i + 1, j)] - p[IX(i - 1, j)]);
                v[IX(i, j)] -= 0.5 * N * (p[IX(i, j + 1)] - p[IX(i, j - 1)]);
            }
        }
        set_bnd(1, u); set_bnd(2, v);
    }
};

// --- MAIN: ВИЗУАЛИЗАЦИЯ ---
int main() {
    // Создаем окно с учетом масштаба
    sf::RenderWindow window(sf::VideoMode((N + 2) * SCALE, (N + 2) * SCALE), "Fluid Simulation 2D");
    window.setFramerateLimit(60);

    GasSolver gas;

    // Текстура, в которую мы будем рисовать пиксели
    sf::Texture texture;
    texture.create(N + 2, N + 2);
    sf::Sprite sprite(texture);
    sprite.setScale(SCALE, SCALE); // Растягиваем маленькую текстуру на все окно

    // Массив пикселей (RGBA - 4 байта на пиксель)
    vector<sf::Uint8> pixels((N + 2) * (N + 2) * 4);

    // Переменные для отслеживания мыши
    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
        }

        // --- 1. ВВОД ПОЛЬЗОВАТЕЛЯ (МЫШЬ) ---
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);
            
            // Переводим координаты экрана в координаты сетки (делим на SCALE)
            int i = mousePos.x / SCALE;
            int j = mousePos.y / SCALE;

            // Ограничение границ
            if (i > 1 && i < N && j > 1 && j < N) {
                // Добавляем газ (плотность)
                gas.addDensity(i, j, 100.0);

                // Вычисляем вектор движения мыши, чтобы "толкнуть" газ
                double forceX = (mousePos.x - prevMousePos.x) * 5.0;
                double forceY = (mousePos.y - prevMousePos.y) * 5.0;
                gas.addVelocity(i, j, forceX, forceY);
            }
            prevMousePos = mousePos;
        } else {
            // Если мышь не нажата, просто обновляем позицию, чтобы не было рывков
            prevMousePos = sf::Mouse::getPosition(window);
        }

        // --- 2. РАСЧЕТ ФИЗИКИ ---
        gas.step();

        // --- 3. РЕНДЕРИНГ ---
        // Преобразуем плотность газа в цвета пикселей
        for (int i = 0; i < N + 2; i++) {
            for (int j = 0; j < N + 2; j++) {
                int index = IX(i, j);
                double d = gas.dens[index];

                // Ограничиваем значение до 255
                int c = (int)d;
                if (c > 255) c = 255;
                if (c < 0) c = 0;

                // Заполняем RGBA (4 байта подряд)
                int pixelIndex = (i + j * (N + 2)) * 4;
                pixels[pixelIndex] = c;     // R
                pixels[pixelIndex + 1] = c; // G
                pixels[pixelIndex + 2] = c; // B
                pixels[pixelIndex + 3] = 255; // Alpha (непрозрачность)
            }
        }

        // Загружаем пиксели в видеокарту
        texture.update(pixels.data());

        window.clear();
        window.draw(sprite);
        window.display();
    }

    return 0;
}