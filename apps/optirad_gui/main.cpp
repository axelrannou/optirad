#include "Application.hpp"
#include "Logger.hpp"

using namespace optirad;

int main(int argc, char* argv[]) {
    Logger::init();
    Logger::info("OptiRad GUI starting...");

    Application app;
    
    if (!app.init()) {
        Logger::error("Failed to initialize application");
        return 1;
    }

    app.run();
    app.shutdown();

    Logger::info("OptiRad GUI finished.");
    return 0;
}
