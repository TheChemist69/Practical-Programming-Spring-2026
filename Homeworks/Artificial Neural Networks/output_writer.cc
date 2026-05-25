#include "output_writer.h"

#include <iomanip>

namespace report {

void write_rule(std::ostream& out) {
  out << "# ------------------------------------------------------------\n";
}

void write_title(std::ostream& out, const char* title) {
  write_rule(out);
  out << "# " << title << "\n";
  write_rule(out);
}

void write_kv(std::ostream& out, const char* key, const char* value) {
  out << "# " << std::left << std::setw(24) << key << ": " << value << "\n";
  out << std::right;
}

void write_kv(std::ostream& out, const char* key, int value) {
  out << "# " << std::left << std::setw(24) << key << ": " << value << "\n";
  out << std::right;
}

void write_kv(std::ostream& out, const char* key, double value) {
  out << "# " << std::left << std::setw(24) << key << ": " << value << "\n";
  out << std::right;
}

void write_interval(std::ostream& out, const char* key, double a, double b) {
  out << "# " << std::left << std::setw(24) << key << ": [" << a << ", " << b
      << "]\n";
  out << std::right;
}

}  // namespace report
