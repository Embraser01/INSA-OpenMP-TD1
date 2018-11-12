#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <unistd.h>
#include <omp.h>

typedef struct Config {
    unsigned int core = 1;
    unsigned int size = 0;
    unsigned int factor = 10;
    bool progressive = false;
} Config;


enum Func {
    addFunc, sumFunc, multFunc
};


typedef std::chrono::high_resolution_clock Clock;


void fillVector(std::vector<double> &v) {
    std::uniform_real_distribution<double> unif(0, 100);
    std::random_device rd;
    std::default_random_engine re(rd());

    auto gen = [&unif, &re]() {
        return unif(re);
    };

    std::generate(v.begin(), v.end(), gen);
}

std::vector<double> generateVector(unsigned int n) {
    std::vector<double> v(n);
    fillVector(v);
    return v;
}

void printVector(std::vector<double> &vector, std::ostream &out, int showLimit = -1) {
    if (showLimit == -1) {
        for (double i : vector) {
            out << i << ", ";
        }
    } else {
        for (int i = 0; i < showLimit; ++i) {
            out << vector[i] << ", ";
        }
    }
}

std::vector<double> add(std::vector<double> &a, std::vector<double> &b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("Vectors must be the same size");
    }
    std::vector<double> sum(a.size());

#pragma omp parallel for shared(sum)
    for (int i = 0; i < a.size(); ++i) {
        sum[i] = a[i] + b[i];
    }

    return sum;
}

double sum(std::vector<double> &v) {
    double acc = 0;

#pragma omp parallel for reduction (*:acc)
    for (int i = 0; i < v.size(); ++i) {
        acc += v[i];
    }

    return acc;
}

Config parseConfig(int argc, char **argv) {
    // Shut GetOpt error messages down (return '?'):
    opterr = 0;
    int opt;

    Config config;

    while ((opt = getopt(argc, argv, "pc:s:f:")) != -1) {
        switch (opt) {
            case 'p':
                config.progressive = true;
                break;
            case 's':
                config.size = static_cast<unsigned int>(std::stoi(optarg));
                break;
            case 'c':
                config.core = static_cast<unsigned int>(std::stoi(optarg));
                break;
            case 'f':
                config.factor = static_cast<unsigned int>(std::stoi(optarg));
                break;
            case '?':
                std::cerr << "Unknown option: '" << char(optopt) << "'!" << std::endl;
                throw 1;
            default:
                break;
        }
    }

    if (config.size == 0) {
        std::cout << "Size must be set with the -s option" << std::endl;
        throw 1;
    }
    if (config.core == 0) {
        std::cout << "Core must be set with the -c option" << std::endl;
        throw 1;
    }

    return config;
}

long computeFun(const Config &config, std::vector<double> &v1, std::vector<double> &v2, const Func &type) {
    long acc = 0;
    for (int i = 0; i < config.factor; ++i) {
//        std::cout << "-- Starting round " << i << std::endl;

        auto start = std::chrono::_V2::system_clock::now();

        switch (type) {
            case addFunc:
                add(v1, v2);
                break;
            case sumFunc:
                sum(v1);
                break;
            case multFunc:
                // TODO Add mult
                break;
        }

        auto end = std::chrono::_V2::system_clock::now();

        acc += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }
    return (acc / config.factor);
}

int main(int argc, char **argv) {
    Config config = parseConfig(argc, argv);

    std::cout << "Size used: " << config.size << std::endl;
    std::cout << "Core used: " << config.core << std::endl;
    std::cout << "Factor used: " << config.factor << std::endl;

    if (config.progressive) {
        // TODO Progressive mode
        //  Increase size with cores and without restarting the app
        Config tempConfig = config;
        tempConfig.core = 1;

        while (tempConfig.core <= config.core) {
            omp_set_num_threads(tempConfig.core);

            std::cout << "\tUsing following config (size, cores): " << tempConfig.size << ", " << tempConfig.core
                      << std::endl;

            std::vector<double> v1 = generateVector(tempConfig.size);
            std::vector<double> v2 = generateVector(tempConfig.size);

            long accAdd = computeFun(tempConfig, v1, v2, addFunc);
            long accSum = computeFun(tempConfig, v1, v2, sumFunc);
            long accMult = computeFun(tempConfig, v1, v2, multFunc);

            std::cout << "Operation took approximately (add, sum, mult) :" << std::endl
                      << accAdd << "ms | " << accSum << "ms | " << accMult << "ms\n";

            // Update config
            tempConfig.size *= 2;
            tempConfig.core *= 2;
        }
    } else {
        omp_set_num_threads(config.core);

        std::vector<double> v1 = generateVector(config.size);
        std::vector<double> v2 = generateVector(config.size);

        long accAdd = computeFun(config, v1, v2, addFunc);
        long accSum = computeFun(config, v1, v2, sumFunc);
        long accMult = computeFun(config, v1, v2, multFunc);
        std::cout << "Operation took approximately (add, sum, mult) :" << std::endl
                  << accAdd << "ms | " << accSum << "ms | " << accMult << "ms\n";
    }

    return 0;
}

