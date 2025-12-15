#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;

// --- НАСТРОЙКИ ВЗРЫВА ---
const int N = 300;           // Хорошее разрешение для Ryzen
const int SCALE = 2;         // Размер окна ~600x600
const double DT = 0.1;       // Шаг времени
const double VORTICITY = 80.0; // ОЧЕНЬ высокая турбулентность для детализации краев

inline int IX(int x, int y) {
    if (x < 0) x = 0; if (x > N + 1) x = N + 1;
    if (y < 0) y = 0; if (y > N + 1) y = N + 1;
    return x + (N + 2) * y;
}

class GasSolver {
public:
    vector<double> u, v, u_prev, v_prev, dens, dens_prev;
    vector<double> curl;

    GasSolver() {
        int size = (N + 2) * (N + 2);
        u.resize(size, 0.0); v.resize(size, 0.0);
        u_prev.resize(size, 0.0); v_prev.resize(size, 0.0);
        dens.resize(size, 0.0); dens_prev.resize(size, 0.0);
        curl.resize(size, 0.0);
    }

    void addDensity(int x, int y, double amount) { dens[IX(x, y)] += amount; }
    void addVelocity(int x, int y, double amountX, double amountY) {
        u[IX(x, y)] += amountX;
        v[IX(x, y)] += amountY;
    }

    void step() {
        vorticityConfinement(u, v);
        
        diffuse(1, u_prev, u, 0.0000); // 0 диффузии для резкости
        diffuse(2, v_prev, v, 0.0000);
        project(u_prev, v_prev, u, v);
        advect(1, u, u_prev, u_prev, v_prev);
        advect(2, v, v_prev, u_prev, v_prev);
        project(u, v, u_prev, v_prev);
        
        diffuse(0, dens_prev, dens, 0.0000);
        advect(0, dens, dens_prev, u, v);

        // Очень медленное остывание (газ расширяется, а не исчезает)
        for(size_t i=0; i<dens.size(); i++) dens[i] *= 0.999;
    }

    // --- ИНИЦИАЛИЗАЦИЯ ВЗРЫВА ---
    void initExplosion() {
        // Очистка
        fill(u.begin(), u.end(), 0.0);
        fill(v.begin(), v.end(), 0.0);
        fill(dens.begin(), dens.end(), 0.0);
        fill(u_prev.begin(), u_prev.end(), 0.0);
        fill(v_prev.begin(), v_prev.end(), 0.0);

        int cx = N / 2;
        int cy = N / 2;
        int radius = N / 10; // Компактное ядро
        double force = 2000.0; // Чудовищная сила

        for (int i = 0; i < N + 2; i++) {
            for (int j = 0; j < N + 2; j++) {
                double dx = i - cx;
                double dy = j - cy;
                double dist = sqrt(dx*dx + dy*dy);

                if (dist < radius) {
                    // 1. Плотность: Очень высокая внутри
                    double d = 1000.0; 
                    addDensity(i, j, d);

                    // 2. Скорость: Радиальная + ШУМ
                    // Шум (rand) критичен для создания "пальцев" неустойчивости
                    double noise = 1.0 + ((rand() % 100) / 50.0 - 1.0) * 0.5; // +/- 50% шума
                    
                    double uVel = (dx / (dist + 0.1)) * force * noise;
                    double vVel = (dy / (dist + 0.1)) * force * noise;

                    addVelocity(i, j, uVel, vVel);
                }
            }
        }
    }

private:
    void vorticityConfinement(vector<double>& u, vector<double>& v) {
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                double du_dy = (u[IX(i, j + 1)] - u[IX(i, j - 1)]) * 0.5;
                double dv_dx = (v[IX(i + 1, j)] - v[IX(i - 1, j)]) * 0.5;
                curl[IX(i, j)] = abs(dv_dx - du_dy);
            }
        }
        for (int i = 1; i <= N; i++) {
            for (int j = 1; j <= N; j++) {
                double dw_dx = (curl[IX(i + 1, j)] - curl[IX(i - 1, j)]) * 0.5;
                double dw_dy = (curl[IX(i, j + 1)] - curl[IX(i, j - 1)]) * 0.5;
                double len = sqrt(dw_dx * dw_dx + dw_dy * dw_dy) + 0.000001;
                u[IX(i, j)] += (dw_dy * -1.0 / len) * curl[IX(i, j)] * VORTICITY * DT;
                v[IX(i, j)] += (dw_dx / len) * curl[IX(i, j)] * VORTICITY * DT;
            }
        }
    }

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
        double a = DT * diff * N * N;
        for (int k = 0; k < 5; k++) {
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
        double dt0 = DT * N;
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
        for (int k = 0; k < 5; k++) {
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

// Цвета взрыва (Гамма-коррекция для яркости)
sf::Color getExplosionColor(double d) {
    if (d <= 1.0) return sf::Color::Black;

    // Сдвигаем диапазон: от 0 до 800
    int val = (int)d;
    if (val > 1000) val = 1000;

    sf::Uint8 r, g, b;

    // 1. Края (Ударная волна) - Темно-красный / Коричневый
    if (val < 100) {
        r = val * 2; g = 0; b = 0;
    }
    // 2. Тело взрыва - Красный -> Оранжевый
    else if (val < 400) {
        r = 255; 
        g = (val - 100) * 0.8; 
        b = 0;
    }
    // 3. Ядро - Желтый -> Белый (Слепящий свет)
    else {
        r = 255;
        g = 255;
        b = (val - 400) * 0.5;
    }

    return sf::Color(r, g, b, 255);
}

int main() {
    sf::RenderWindow window(sf::VideoMode((N + 2) * SCALE, (N + 2) * SCALE), "Shockwave Simulation");
    window.setFramerateLimit(60);

    GasSolver gas;
    
    // БУМ! Запускаем сразу при старте
    gas.initExplosion();

    sf::Texture texture;
    texture.create(N + 2, N + 2);
    sf::Sprite sprite(texture);
    sprite.setScale(SCALE, SCALE);
    vector<sf::Uint8> pixels((N + 2) * (N + 2) * 4);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) window.close();
            // Перезапуск взрыва на Пробел
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Space) {
                gas.initExplosion();
            }
        }

        gas.step();

        for (int i = 0; i < N + 2; i++) {
            for (int j = 0; j < N + 2; j++) {
                int idx = IX(i, j);
                sf::Color c = getExplosionColor(gas.dens[idx]);
                int pIdx = idx * 4;
                pixels[pIdx] = c.r;
                pixels[pIdx+1] = c.g;
                pixels[pIdx+2] = c.b;
                pixels[pIdx+3] = 255;
            }
        }

        texture.update(pixels.data());
        window.clear();
        window.draw(sprite);
        window.display();
    }
    return 0;
}