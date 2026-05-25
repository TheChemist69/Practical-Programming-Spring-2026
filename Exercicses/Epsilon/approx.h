#pragma once

// Returns true if a and b are approximately equal (absolute or relative tolerance)
bool approx(double a, double b, double acc = 1e-9, double eps = 1e-9);
