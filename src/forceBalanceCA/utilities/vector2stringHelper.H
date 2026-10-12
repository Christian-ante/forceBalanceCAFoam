// HELPER FUNCTION TO PRINT VECTOR TO TXT FILE
#include <sstream>
std::string vectorToString(const Foam::vector& vec)
{
    std::ostringstream oss;
    oss << "(" << vec.x() << ", " << vec.y() << ", " << vec.z() << ")";
    return oss.str();
}
