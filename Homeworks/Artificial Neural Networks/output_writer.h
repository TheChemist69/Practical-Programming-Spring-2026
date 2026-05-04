#ifndef OUTPUT_WRITER_H
#define OUTPUT_WRITER_H

#include <ostream>

namespace report {

void write_rule(std::ostream& out);
void write_title(std::ostream& out, const char* title);
void write_kv(std::ostream& out, const char* key, const char* value);
void write_kv(std::ostream& out, const char* key, int value);
void write_kv(std::ostream& out, const char* key, double value);
void write_interval(std::ostream& out, const char* key, double a, double b);

}  // namespace report

#endif  // OUTPUT_WRITER_H
