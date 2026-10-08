#include "validation/DoseRawIO.hpp"

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace optirad {

namespace {
constexpr char kMagic[8] = {'O', 'D', 'O', 'S', 'E', '1', 0, 0};
}

void DoseRawIO::write(const std::string& path, const DoseMatrix& dose) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("DoseRawIO: cannot open for writing: " + path);

    const auto dims = dose.getGrid().getDimensions();
    const int32_t n[3] = {int32_t(dims[0]), int32_t(dims[1]), int32_t(dims[2])};
    const Vec3 sp = dose.getGrid().getSpacing();
    const Vec3 org = dose.getGrid().getOrigin();

    f.write(kMagic, 8);
    f.write(reinterpret_cast<const char*>(n), sizeof(n));
    f.write(reinterpret_cast<const char*>(sp.data()), sizeof(double) * 3);
    f.write(reinterpret_cast<const char*>(org.data()), sizeof(double) * 3);
    f.write(reinterpret_cast<const char*>(dose.data()), std::streamsize(sizeof(double) * dose.size()));
    if (!f) throw std::runtime_error("DoseRawIO: write failed: " + path);
}

std::shared_ptr<DoseMatrix> DoseRawIO::read(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("DoseRawIO: cannot open: " + path);

    char magic[8];
    int32_t n[3];
    double sp[3], org[3];
    f.read(magic, 8);
    if (!f || std::memcmp(magic, kMagic, 8) != 0)
        throw std::runtime_error("DoseRawIO: not an ODOSE1 file: " + path);
    f.read(reinterpret_cast<char*>(n), sizeof(n));
    f.read(reinterpret_cast<char*>(sp), sizeof(sp));
    f.read(reinterpret_cast<char*>(org), sizeof(org));
    if (!f || n[0] <= 0 || n[1] <= 0 || n[2] <= 0)
        throw std::runtime_error("DoseRawIO: invalid header: " + path);

    Grid g;
    g.setDimensions(size_t(n[0]), size_t(n[1]), size_t(n[2]));
    g.setSpacing(sp[0], sp[1], sp[2]);
    g.setOrigin(Vec3{org[0], org[1], org[2]});

    auto dose = std::make_shared<DoseMatrix>();
    dose->setGrid(g);
    dose->allocate();
    f.read(reinterpret_cast<char*>(dose->data()), std::streamsize(sizeof(double) * dose->size()));
    if (!f) throw std::runtime_error("DoseRawIO: truncated data: " + path);
    return dose;
}

} // namespace optirad
