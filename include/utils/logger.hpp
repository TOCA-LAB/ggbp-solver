#pragma once
#include <iostream>
#include <sstream>

#include "utils/fix_c++_print.hpp"

class Logger {
public:
    int verbosity_level = 0;

    void set_level(int v) { verbosity_level = v; }

    class Entry {
    public:
        Entry(int level, bool enabled) 
            : level_(level), enabled_(enabled) {}

        // operador << para qualquer tipo
        template <typename T>
        Entry& operator<<(const T& value) {
            if (enabled_) {
                stream_ << value;
            }
            return *this;
        }

        // manipuladores como std::endl
        Entry& operator<<(std::ostream& (*manip)(std::ostream&)) {
            if (enabled_) {
                manip(stream_);
            }
            return *this;
        }

        ~Entry() {
            if (enabled_) {
                std::cerr << "[INFO]\t";
                std::cerr << stream_.str() << "\n";
                std::cerr.flush();
            }
        }

    private:
        int level_;
        bool enabled_;
        std::ostringstream stream_;
    };

    Entry info(int level) const {
        return Entry(level, verbosity_level >= level);
    }
};

inline Logger logger;
