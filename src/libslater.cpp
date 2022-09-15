#include "logger.h"

namespace slater {

    void libslater_global_init()
    {
        logger()->info("Slater library initialized");
    }


    void libslater_global_cleanup()
    {
        logger()->info("Slater library cleaned up");
    }

}