// Exercise the production gravity model on a deterministic density field.
#include "ReadParameters.h"
#include "Cosmology_LCDM.h"
#include "Cosmology_HiCOLA.h"
#include "GravityModel_GR.h"
#include "GravityModel_HiCOLA.h"
#include "GravityModel_kmouflage.h"
#include <fstream>
#include <iomanip>

int main(int argc, char **argv) {
    if (argc != 3) throw std::runtime_error("Usage: force_probe run.lua output-prefix");
    ParameterMap param;
    read_parameterfile(param, argv[1]);
    std::shared_ptr<Cosmology> cosmo;
    if (param.get<std::string>("cosmology_model") == "LCDM")
        cosmo = std::make_shared<CosmologyLCDM>();
    else
        cosmo = std::make_shared<CosmologyHiCOLA>();
    cosmo->read_parameters(param);
    if (cosmo->get_name() == "HiCOLA") cosmo->read_and_spline_expansion();
    cosmo->init();
    std::shared_ptr<GravityModel<3>> grav;
    auto name = param.get<std::string>("gravity_model");
    if (name == "GR") grav = std::make_shared<GravityModelGR<3>>(cosmo);
    else if (name == "HiCOLA") grav = std::make_shared<GravityModelHiCOLA<3>>(cosmo);
    else grav = std::make_shared<GravityModelkmouflage<3>>(cosmo);
    grav->read_parameters(param);
    if (name != "GR") grav->read_and_spline_chi_over_delta_chi();
    grav->init();
    FML::GRID::FFTWGrid<3> density(16, 0, 0);
    for (auto index : density.get_real_range()) {
        auto pos = density.get_real_position(density.get_coord_from_index(index));
        density.set_real_from_index(index, .2 * (std::sin(2*M_PI*pos[0]) +
            std::cos(4*M_PI*pos[1]) + std::sin(6*M_PI*pos[2])));
    }
    density.fftw_r2c();
    std::array<FML::GRID::FFTWGrid<3>,3> force;
    grav->compute_force(.8, param.get<double>("simulation_boxsize") * grav->H0_hmpc, density, "CIC", force);
    std::ofstream out(std::string(argv[2]) + "." + std::to_string(FML::ThisTask) + ".txt");
    out << std::setprecision(17);
    for (auto index : force[0].get_real_range()) {
        auto pos = force[0].get_real_position(force[0].get_coord_from_index(index));
        for (auto x : pos) out << x << ' ';
        for (auto &component : force) out << component.get_real_from_index(index) << ' ';
        out << '\n';
    }
}
