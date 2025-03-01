// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_combintersection.hpp"
#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_generated_6890()
{
  Cut::MeshIntersection intersection;
  intersection.get_options().init_for_cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  std::vector<double> lsvs(8);
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.141497;
    tri3_xyze(1, 0) = 0.0379141;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(8335);
    tri3_xyze(0, 1) = 0.146489;
    tri3_xyze(1, 1) = -5.78249e-17;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = 0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-673);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0996251;
    tri3_xyze(1, 0) = -0.0266945;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5317);
    tri3_xyze(0, 1) = 0.105642;
    tri3_xyze(1, 1) = -0.0283068;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 2) = 0.0973263;
    tri3_xyze(1, 2) = -0.0403139;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-565);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105642;
    tri3_xyze(1, 0) = -0.0283068;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 1) = 0.0947164;
    tri3_xyze(1, 1) = -0.0546845;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5069);
    tri3_xyze(0, 2) = 0.0973263;
    tri3_xyze(1, 2) = -0.0403139;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-565);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105642;
    tri3_xyze(1, 0) = -0.0283068;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 1) = 0.111517;
    tri3_xyze(1, 1) = -0.0298809;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 2) = 0.102965;
    tri3_xyze(1, 2) = -0.0426494;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-566);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.111517;
    tri3_xyze(1, 0) = -0.0298809;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 1) = 0.0999834;
    tri3_xyze(1, 1) = -0.0577254;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5071);
    tri3_xyze(0, 2) = 0.102965;
    tri3_xyze(1, 2) = -0.0426494;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-566);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0999834;
    tri3_xyze(1, 0) = -0.0577254;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5071);
    tri3_xyze(0, 1) = 0.0947164;
    tri3_xyze(1, 1) = -0.0546845;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5069);
    tri3_xyze(0, 2) = 0.102965;
    tri3_xyze(1, 2) = -0.0426494;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-566);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0947164;
    tri3_xyze(1, 0) = -0.0546845;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5069);
    tri3_xyze(0, 1) = 0.105642;
    tri3_xyze(1, 1) = -0.0283068;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 2) = 0.102965;
    tri3_xyze(1, 2) = -0.0426494;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-566);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.111517;
    tri3_xyze(1, 0) = -0.0298809;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 1) = 0.117156;
    tri3_xyze(1, 1) = -0.0313919;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 2) = 0.108424;
    tri3_xyze(1, 2) = -0.0449107;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-567);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.117156;
    tri3_xyze(1, 0) = -0.0313919;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 1) = 0.105039;
    tri3_xyze(1, 1) = -0.0606445;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5073);
    tri3_xyze(0, 2) = 0.108424;
    tri3_xyze(1, 2) = -0.0449107;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-567);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105039;
    tri3_xyze(1, 0) = -0.0606445;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5073);
    tri3_xyze(0, 1) = 0.0999834;
    tri3_xyze(1, 1) = -0.0577254;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5071);
    tri3_xyze(0, 2) = 0.108424;
    tri3_xyze(1, 2) = -0.0449107;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-567);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0999834;
    tri3_xyze(1, 0) = -0.0577254;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5071);
    tri3_xyze(0, 1) = 0.111517;
    tri3_xyze(1, 1) = -0.0298809;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 2) = 0.108424;
    tri3_xyze(1, 2) = -0.0449107;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-567);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.117156;
    tri3_xyze(1, 0) = -0.0313919;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 1) = 0.122471;
    tri3_xyze(1, 1) = -0.032816;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 2) = 0.113618;
    tri3_xyze(1, 2) = -0.047062;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-568);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.122471;
    tri3_xyze(1, 0) = -0.032816;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 1) = 0.109805;
    tri3_xyze(1, 1) = -0.0633957;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5075);
    tri3_xyze(0, 2) = 0.113618;
    tri3_xyze(1, 2) = -0.047062;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-568);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.109805;
    tri3_xyze(1, 0) = -0.0633957;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5075);
    tri3_xyze(0, 1) = 0.105039;
    tri3_xyze(1, 1) = -0.0606445;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5073);
    tri3_xyze(0, 2) = 0.113618;
    tri3_xyze(1, 2) = -0.047062;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-568);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105039;
    tri3_xyze(1, 0) = -0.0606445;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5073);
    tri3_xyze(0, 1) = 0.117156;
    tri3_xyze(1, 1) = -0.0313919;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 2) = 0.113618;
    tri3_xyze(1, 2) = -0.047062;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-568);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.122471;
    tri3_xyze(1, 0) = -0.032816;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 1) = 0.127378;
    tri3_xyze(1, 1) = -0.0341308;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 2) = 0.118464;
    tri3_xyze(1, 2) = -0.0490695;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-569);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.127378;
    tri3_xyze(1, 0) = -0.0341308;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 1) = 0.114204;
    tri3_xyze(1, 1) = -0.0659356;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5077);
    tri3_xyze(0, 2) = 0.118464;
    tri3_xyze(1, 2) = -0.0490695;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-569);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.114204;
    tri3_xyze(1, 0) = -0.0659356;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5077);
    tri3_xyze(0, 1) = 0.109805;
    tri3_xyze(1, 1) = -0.0633957;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5075);
    tri3_xyze(0, 2) = 0.118464;
    tri3_xyze(1, 2) = -0.0490695;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-569);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.109805;
    tri3_xyze(1, 0) = -0.0633957;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5075);
    tri3_xyze(0, 1) = 0.122471;
    tri3_xyze(1, 1) = -0.032816;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 2) = 0.118464;
    tri3_xyze(1, 2) = -0.0490695;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-569);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.127378;
    tri3_xyze(1, 0) = -0.0341308;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 1) = 0.131799;
    tri3_xyze(1, 1) = -0.0353155;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 2) = 0.122887;
    tri3_xyze(1, 2) = -0.0509015;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-570);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131799;
    tri3_xyze(1, 0) = -0.0353155;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 1) = 0.118168;
    tri3_xyze(1, 1) = -0.0682242;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5079);
    tri3_xyze(0, 2) = 0.122887;
    tri3_xyze(1, 2) = -0.0509015;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-570);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.114204;
    tri3_xyze(1, 0) = -0.0659356;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5077);
    tri3_xyze(0, 1) = 0.127378;
    tri3_xyze(1, 1) = -0.0341308;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 2) = 0.122887;
    tri3_xyze(1, 2) = -0.0509015;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-570);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131799;
    tri3_xyze(1, 0) = -0.0353155;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 1) = 0.135665;
    tri3_xyze(1, 1) = -0.0363514;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 2) = 0.126816;
    tri3_xyze(1, 2) = -0.0525291;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-571);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = -0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 1) = 0.121634;
    tri3_xyze(1, 1) = -0.0702254;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5081);
    tri3_xyze(0, 2) = 0.126816;
    tri3_xyze(1, 2) = -0.0525291;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-571);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.118168;
    tri3_xyze(1, 0) = -0.0682242;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5079);
    tri3_xyze(0, 1) = 0.131799;
    tri3_xyze(1, 1) = -0.0353155;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 2) = 0.126816;
    tri3_xyze(1, 2) = -0.0525291;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-571);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = -0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 1) = 0.138915;
    tri3_xyze(1, 1) = -0.0372221;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 2) = 0.13019;
    tri3_xyze(1, 2) = -0.0539266;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-572);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.138915;
    tri3_xyze(1, 0) = -0.0372221;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 1) = 0.124548;
    tri3_xyze(1, 1) = -0.0719077;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5083);
    tri3_xyze(0, 2) = 0.13019;
    tri3_xyze(1, 2) = -0.0539266;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-572);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121634;
    tri3_xyze(1, 0) = -0.0702254;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5081);
    tri3_xyze(0, 1) = 0.135665;
    tri3_xyze(1, 1) = -0.0363514;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 2) = 0.13019;
    tri3_xyze(1, 2) = -0.0539266;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-572);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.138915;
    tri3_xyze(1, 0) = -0.0372221;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 1) = 0.141497;
    tri3_xyze(1, 1) = -0.0379141;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 2) = 0.132956;
    tri3_xyze(1, 2) = -0.0550721;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-573);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.141497;
    tri3_xyze(1, 0) = -0.0379141;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 1) = 0.126863;
    tri3_xyze(1, 1) = -0.0732444;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5085);
    tri3_xyze(0, 2) = 0.132956;
    tri3_xyze(1, 2) = -0.0550721;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-573);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.124548;
    tri3_xyze(1, 0) = -0.0719077;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5083);
    tri3_xyze(0, 1) = 0.138915;
    tri3_xyze(1, 1) = -0.0372221;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 2) = 0.132956;
    tri3_xyze(1, 2) = -0.0550721;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-573);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.141497;
    tri3_xyze(1, 0) = -0.0379141;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 1) = 0.143372;
    tri3_xyze(1, 1) = -0.0384163;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 2) = 0.135069;
    tri3_xyze(1, 2) = -0.0559473;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-574);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143372;
    tri3_xyze(1, 0) = -0.0384163;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 1) = 0.128543;
    tri3_xyze(1, 1) = -0.0742146;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5087);
    tri3_xyze(0, 2) = 0.135069;
    tri3_xyze(1, 2) = -0.0559473;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-574);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.126863;
    tri3_xyze(1, 0) = -0.0732444;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5085);
    tri3_xyze(0, 1) = 0.141497;
    tri3_xyze(1, 1) = -0.0379141;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 2) = 0.135069;
    tri3_xyze(1, 2) = -0.0559473;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-574);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143372;
    tri3_xyze(1, 0) = -0.0384163;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 1) = 0.144508;
    tri3_xyze(1, 1) = -0.0387208;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 2) = 0.136496;
    tri3_xyze(1, 2) = -0.0565386;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-575);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144508;
    tri3_xyze(1, 0) = -0.0387208;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 1) = 0.129562;
    tri3_xyze(1, 1) = -0.0748029;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5089);
    tri3_xyze(0, 2) = 0.136496;
    tri3_xyze(1, 2) = -0.0565386;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-575);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.128543;
    tri3_xyze(1, 0) = -0.0742146;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5087);
    tri3_xyze(0, 1) = 0.143372;
    tri3_xyze(1, 1) = -0.0384163;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 2) = 0.136496;
    tri3_xyze(1, 2) = -0.0565386;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-575);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144508;
    tri3_xyze(1, 0) = -0.0387208;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = -0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 2) = 0.137216;
    tri3_xyze(1, 2) = -0.0568366;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-576);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144889;
    tri3_xyze(1, 0) = -0.0388229;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 1) = 0.129904;
    tri3_xyze(1, 1) = -0.075;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5091);
    tri3_xyze(0, 2) = 0.137216;
    tri3_xyze(1, 2) = -0.0568366;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-576);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.129562;
    tri3_xyze(1, 0) = -0.0748029;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5089);
    tri3_xyze(0, 1) = 0.144508;
    tri3_xyze(1, 1) = -0.0387208;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 2) = 0.137216;
    tri3_xyze(1, 2) = -0.0568366;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-576);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144889;
    tri3_xyze(1, 0) = -0.0388229;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 1) = 0.144508;
    tri3_xyze(1, 1) = -0.0387208;
    tri3_xyze(2, 1) = 0.806267;
    nids.push_back(5343);
    tri3_xyze(0, 2) = 0.137216;
    tri3_xyze(1, 2) = -0.0568366;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-577);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.129904;
    tri3_xyze(1, 0) = -0.075;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5091);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = -0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 2) = 0.137216;
    tri3_xyze(1, 2) = -0.0568366;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-577);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0968605;
    tri3_xyze(1, 0) = -1.5677e-16;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5565);
    tri3_xyze(0, 1) = 0.10314;
    tri3_xyze(1, 1) = -1.5677e-16;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 2) = 0.0982963;
    tri3_xyze(1, 2) = -0.012941;
    tri3_xyze(2, 2) = 0.750099;
    nids.push_back(-614);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10314;
    tri3_xyze(1, 0) = -1.5677e-16;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 1) = 0.0996251;
    tri3_xyze(1, 1) = -0.0266945;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5317);
    tri3_xyze(0, 2) = 0.0982963;
    tri3_xyze(1, 2) = -0.012941;
    tri3_xyze(2, 2) = 0.750099;
    nids.push_back(-614);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10314;
    tri3_xyze(1, 0) = -1.5677e-16;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 1) = 0.109369;
    tri3_xyze(1, 1) = -1.54297e-16;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = -0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-615);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.109369;
    tri3_xyze(1, 0) = -1.54297e-16;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 1) = 0.105642;
    tri3_xyze(1, 1) = -0.0283068;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = -0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-615);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105642;
    tri3_xyze(1, 0) = -0.0283068;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 1) = 0.0996251;
    tri3_xyze(1, 1) = -0.0266945;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5317);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = -0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-615);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0996251;
    tri3_xyze(1, 0) = -0.0266945;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5317);
    tri3_xyze(0, 1) = 0.10314;
    tri3_xyze(1, 1) = -1.5677e-16;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = -0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-615);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.109369;
    tri3_xyze(1, 0) = -1.54297e-16;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 1) = 0.115451;
    tri3_xyze(1, 1) = -1.49392e-16;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = -0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-616);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.115451;
    tri3_xyze(1, 0) = -1.49392e-16;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 1) = 0.111517;
    tri3_xyze(1, 1) = -0.0298809;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = -0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-616);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.111517;
    tri3_xyze(1, 0) = -0.0298809;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 1) = 0.105642;
    tri3_xyze(1, 1) = -0.0283068;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = -0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-616);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105642;
    tri3_xyze(1, 0) = -0.0283068;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5319);
    tri3_xyze(0, 1) = 0.109369;
    tri3_xyze(1, 1) = -1.54297e-16;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = -0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-616);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.115451;
    tri3_xyze(1, 0) = -1.49392e-16;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 1) = 0.121289;
    tri3_xyze(1, 1) = -1.4213e-16;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = -0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-617);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121289;
    tri3_xyze(1, 0) = -1.4213e-16;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 1) = 0.117156;
    tri3_xyze(1, 1) = -0.0313919;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = -0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-617);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.117156;
    tri3_xyze(1, 0) = -0.0313919;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 1) = 0.111517;
    tri3_xyze(1, 1) = -0.0298809;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = -0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-617);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.111517;
    tri3_xyze(1, 0) = -0.0298809;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5321);
    tri3_xyze(0, 1) = 0.115451;
    tri3_xyze(1, 1) = -1.49392e-16;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = -0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-617);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121289;
    tri3_xyze(1, 0) = -1.4213e-16;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 1) = 0.126791;
    tri3_xyze(1, 1) = -1.32627e-16;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = -0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-618);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.126791;
    tri3_xyze(1, 0) = -1.32627e-16;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 1) = 0.122471;
    tri3_xyze(1, 1) = -0.032816;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = -0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-618);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.122471;
    tri3_xyze(1, 0) = -0.032816;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 1) = 0.117156;
    tri3_xyze(1, 1) = -0.0313919;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = -0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-618);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.117156;
    tri3_xyze(1, 0) = -0.0313919;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5323);
    tri3_xyze(0, 1) = 0.121289;
    tri3_xyze(1, 1) = -1.4213e-16;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = -0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-618);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.126791;
    tri3_xyze(1, 0) = -1.32627e-16;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 1) = 0.131871;
    tri3_xyze(1, 1) = -1.21032e-16;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = -0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-619);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131871;
    tri3_xyze(1, 0) = -1.21032e-16;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 1) = 0.127378;
    tri3_xyze(1, 1) = -0.0341308;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = -0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-619);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.127378;
    tri3_xyze(1, 0) = -0.0341308;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 1) = 0.122471;
    tri3_xyze(1, 1) = -0.032816;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = -0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-619);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.122471;
    tri3_xyze(1, 0) = -0.032816;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5325);
    tri3_xyze(0, 1) = 0.126791;
    tri3_xyze(1, 1) = -1.32627e-16;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = -0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-619);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131871;
    tri3_xyze(1, 0) = -1.21032e-16;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 1) = 0.136448;
    tri3_xyze(1, 1) = -1.07528e-16;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = -0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-620);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.136448;
    tri3_xyze(1, 0) = -1.07528e-16;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 1) = 0.131799;
    tri3_xyze(1, 1) = -0.0353155;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = -0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-620);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131799;
    tri3_xyze(1, 0) = -0.0353155;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 1) = 0.127378;
    tri3_xyze(1, 1) = -0.0341308;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = -0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-620);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.127378;
    tri3_xyze(1, 0) = -0.0341308;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5327);
    tri3_xyze(0, 1) = 0.131871;
    tri3_xyze(1, 1) = -1.21032e-16;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = -0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-620);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.136448;
    tri3_xyze(1, 0) = -1.07528e-16;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 1) = 0.140451;
    tri3_xyze(1, 1) = -9.23291e-17;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = -0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-621);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.140451;
    tri3_xyze(1, 0) = -9.23291e-17;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 1) = 0.135665;
    tri3_xyze(1, 1) = -0.0363514;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = -0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-621);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = -0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 1) = 0.131799;
    tri3_xyze(1, 1) = -0.0353155;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = -0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-621);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131799;
    tri3_xyze(1, 0) = -0.0353155;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5329);
    tri3_xyze(0, 1) = 0.136448;
    tri3_xyze(1, 1) = -1.07528e-16;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = -0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-621);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.140451;
    tri3_xyze(1, 0) = -9.23291e-17;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 1) = 0.143815;
    tri3_xyze(1, 1) = -7.56737e-17;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = -0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-622);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143815;
    tri3_xyze(1, 0) = -7.56737e-17;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 1) = 0.138915;
    tri3_xyze(1, 1) = -0.0372221;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = -0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-622);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.138915;
    tri3_xyze(1, 0) = -0.0372221;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 1) = 0.135665;
    tri3_xyze(1, 1) = -0.0363514;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = -0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-622);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = -0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5331);
    tri3_xyze(0, 1) = 0.140451;
    tri3_xyze(1, 1) = -9.23291e-17;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = -0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-622);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143815;
    tri3_xyze(1, 0) = -7.56737e-17;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 1) = 0.146489;
    tri3_xyze(1, 1) = -5.78249e-17;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = -0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-623);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.146489;
    tri3_xyze(1, 0) = -5.78249e-17;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 1) = 0.141497;
    tri3_xyze(1, 1) = -0.0379141;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = -0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-623);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.141497;
    tri3_xyze(1, 0) = -0.0379141;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 1) = 0.138915;
    tri3_xyze(1, 1) = -0.0372221;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = -0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-623);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.138915;
    tri3_xyze(1, 0) = -0.0372221;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5333);
    tri3_xyze(0, 1) = 0.143815;
    tri3_xyze(1, 1) = -7.56737e-17;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = -0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-623);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.146489;
    tri3_xyze(1, 0) = -5.78249e-17;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 1) = 0.148429;
    tri3_xyze(1, 1) = -3.90641e-17;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = -0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-624);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.148429;
    tri3_xyze(1, 0) = -3.90641e-17;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 1) = 0.143372;
    tri3_xyze(1, 1) = -0.0384163;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = -0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-624);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143372;
    tri3_xyze(1, 0) = -0.0384163;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 1) = 0.141497;
    tri3_xyze(1, 1) = -0.0379141;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = -0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-624);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.141497;
    tri3_xyze(1, 0) = -0.0379141;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5335);
    tri3_xyze(0, 1) = 0.146489;
    tri3_xyze(1, 1) = -5.78249e-17;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = -0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-624);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.148429;
    tri3_xyze(1, 0) = -3.90641e-17;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 1) = 0.149606;
    tri3_xyze(1, 1) = -1.96873e-17;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = -0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-625);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.149606;
    tri3_xyze(1, 0) = -1.96873e-17;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 1) = 0.144508;
    tri3_xyze(1, 1) = -0.0387208;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = -0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-625);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144508;
    tri3_xyze(1, 0) = -0.0387208;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 1) = 0.143372;
    tri3_xyze(1, 1) = -0.0384163;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = -0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-625);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143372;
    tri3_xyze(1, 0) = -0.0384163;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5337);
    tri3_xyze(0, 1) = 0.148429;
    tri3_xyze(1, 1) = -3.90641e-17;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = -0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-625);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.149606;
    tri3_xyze(1, 0) = -1.96873e-17;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 1) = 0.15;
    tri3_xyze(1, 1) = 0;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-626);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.15;
    tri3_xyze(1, 0) = 0;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = -0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-626);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144889;
    tri3_xyze(1, 0) = -0.0388229;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 1) = 0.144508;
    tri3_xyze(1, 1) = -0.0387208;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-626);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144508;
    tri3_xyze(1, 0) = -0.0387208;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5339);
    tri3_xyze(0, 1) = 0.149606;
    tri3_xyze(1, 1) = -1.96873e-17;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-626);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.15;
    tri3_xyze(1, 0) = 0;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 1) = 0.149606;
    tri3_xyze(1, 1) = 1.96873e-17;
    tri3_xyze(2, 1) = 0.806267;
    nids.push_back(5593);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-627);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144508;
    tri3_xyze(1, 0) = -0.0387208;
    tri3_xyze(2, 0) = 0.806267;
    nids.push_back(5343);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = -0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-627);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144889;
    tri3_xyze(1, 0) = -0.0388229;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5341);
    tri3_xyze(0, 1) = 0.15;
    tri3_xyze(1, 1) = 0;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = -0.0193859;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-627);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.0996251;
    tri3_xyze(1, 0) = 0.0266945;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(8317);
    tri3_xyze(0, 1) = 0.10314;
    tri3_xyze(1, 1) = -1.5677e-16;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 2) = 0.0982963;
    tri3_xyze(1, 2) = 0.012941;
    tri3_xyze(2, 2) = 0.750099;
    nids.push_back(-664);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10314;
    tri3_xyze(1, 0) = -1.5677e-16;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 1) = 0.0968605;
    tri3_xyze(1, 1) = -1.5677e-16;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5565);
    tri3_xyze(0, 2) = 0.0982963;
    tri3_xyze(1, 2) = 0.012941;
    tri3_xyze(2, 2) = 0.750099;
    nids.push_back(-664);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.105642;
    tri3_xyze(1, 0) = 0.0283068;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(8319);
    tri3_xyze(0, 1) = 0.109369;
    tri3_xyze(1, 1) = -1.54297e-16;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = 0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-665);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.109369;
    tri3_xyze(1, 0) = -1.54297e-16;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 1) = 0.10314;
    tri3_xyze(1, 1) = -1.5677e-16;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = 0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-665);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.10314;
    tri3_xyze(1, 0) = -1.5677e-16;
    tri3_xyze(2, 0) = 0.750099;
    nids.push_back(5567);
    tri3_xyze(0, 1) = 0.0996251;
    tri3_xyze(1, 1) = 0.0266945;
    tri3_xyze(2, 1) = 0.750099;
    nids.push_back(8317);
    tri3_xyze(0, 2) = 0.104444;
    tri3_xyze(1, 2) = 0.0137503;
    tri3_xyze(2, 2) = 0.750492;
    nids.push_back(-665);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.111517;
    tri3_xyze(1, 0) = 0.0298809;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(8321);
    tri3_xyze(0, 1) = 0.115451;
    tri3_xyze(1, 1) = -1.49392e-16;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = 0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-666);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.115451;
    tri3_xyze(1, 0) = -1.49392e-16;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 1) = 0.109369;
    tri3_xyze(1, 1) = -1.54297e-16;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = 0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-666);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.109369;
    tri3_xyze(1, 0) = -1.54297e-16;
    tri3_xyze(2, 0) = 0.750886;
    nids.push_back(5569);
    tri3_xyze(0, 1) = 0.105642;
    tri3_xyze(1, 1) = 0.0283068;
    tri3_xyze(2, 1) = 0.750886;
    nids.push_back(8319);
    tri3_xyze(0, 2) = 0.110495;
    tri3_xyze(1, 2) = 0.0145469;
    tri3_xyze(2, 2) = 0.751666;
    nids.push_back(-666);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.117156;
    tri3_xyze(1, 0) = 0.0313919;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(8323);
    tri3_xyze(0, 1) = 0.121289;
    tri3_xyze(1, 1) = -1.4213e-16;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = 0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-667);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121289;
    tri3_xyze(1, 0) = -1.4213e-16;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 1) = 0.115451;
    tri3_xyze(1, 1) = -1.49392e-16;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = 0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-667);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.115451;
    tri3_xyze(1, 0) = -1.49392e-16;
    tri3_xyze(2, 0) = 0.752447;
    nids.push_back(5571);
    tri3_xyze(0, 1) = 0.111517;
    tri3_xyze(1, 1) = 0.0298809;
    tri3_xyze(2, 1) = 0.752447;
    nids.push_back(8321);
    tri3_xyze(0, 2) = 0.116353;
    tri3_xyze(1, 2) = 0.0153182;
    tri3_xyze(2, 2) = 0.753603;
    nids.push_back(-667);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.122471;
    tri3_xyze(1, 0) = 0.032816;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(8325);
    tri3_xyze(0, 1) = 0.126791;
    tri3_xyze(1, 1) = -1.32627e-16;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = 0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-668);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.126791;
    tri3_xyze(1, 0) = -1.32627e-16;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 1) = 0.121289;
    tri3_xyze(1, 1) = -1.4213e-16;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = 0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-668);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.121289;
    tri3_xyze(1, 0) = -1.4213e-16;
    tri3_xyze(2, 0) = 0.754759;
    nids.push_back(5573);
    tri3_xyze(0, 1) = 0.117156;
    tri3_xyze(1, 1) = 0.0313919;
    tri3_xyze(2, 1) = 0.754759;
    nids.push_back(8323);
    tri3_xyze(0, 2) = 0.121927;
    tri3_xyze(1, 2) = 0.016052;
    tri3_xyze(2, 2) = 0.756271;
    nids.push_back(-668);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.127378;
    tri3_xyze(1, 0) = 0.0341308;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(8327);
    tri3_xyze(0, 1) = 0.131871;
    tri3_xyze(1, 1) = -1.21032e-16;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = 0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-669);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131871;
    tri3_xyze(1, 0) = -1.21032e-16;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 1) = 0.126791;
    tri3_xyze(1, 1) = -1.32627e-16;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = 0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-669);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.126791;
    tri3_xyze(1, 0) = -1.32627e-16;
    tri3_xyze(2, 0) = 0.757784;
    nids.push_back(5575);
    tri3_xyze(0, 1) = 0.122471;
    tri3_xyze(1, 1) = 0.032816;
    tri3_xyze(2, 1) = 0.757784;
    nids.push_back(8325);
    tri3_xyze(0, 2) = 0.127128;
    tri3_xyze(1, 2) = 0.0167367;
    tri3_xyze(2, 2) = 0.759629;
    nids.push_back(-669);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131799;
    tri3_xyze(1, 0) = 0.0353155;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(8329);
    tri3_xyze(0, 1) = 0.136448;
    tri3_xyze(1, 1) = -1.07528e-16;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = 0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-670);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.136448;
    tri3_xyze(1, 0) = -1.07528e-16;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 1) = 0.131871;
    tri3_xyze(1, 1) = -1.21032e-16;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = 0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-670);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.131871;
    tri3_xyze(1, 0) = -1.21032e-16;
    tri3_xyze(2, 0) = 0.761474;
    nids.push_back(5577);
    tri3_xyze(0, 1) = 0.127378;
    tri3_xyze(1, 1) = 0.0341308;
    tri3_xyze(2, 1) = 0.761474;
    nids.push_back(8327);
    tri3_xyze(0, 2) = 0.131874;
    tri3_xyze(1, 2) = 0.0173616;
    tri3_xyze(2, 2) = 0.763623;
    nids.push_back(-670);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.135665;
    tri3_xyze(1, 0) = 0.0363514;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(8331);
    tri3_xyze(0, 1) = 0.140451;
    tri3_xyze(1, 1) = -9.23291e-17;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = 0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-671);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.140451;
    tri3_xyze(1, 0) = -9.23291e-17;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 1) = 0.136448;
    tri3_xyze(1, 1) = -1.07528e-16;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = 0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-671);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.136448;
    tri3_xyze(1, 0) = -1.07528e-16;
    tri3_xyze(2, 0) = 0.765773;
    nids.push_back(5579);
    tri3_xyze(0, 1) = 0.131799;
    tri3_xyze(1, 1) = 0.0353155;
    tri3_xyze(2, 1) = 0.765773;
    nids.push_back(8329);
    tri3_xyze(0, 2) = 0.136091;
    tri3_xyze(1, 2) = 0.0179167;
    tri3_xyze(2, 2) = 0.768192;
    nids.push_back(-671);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.138915;
    tri3_xyze(1, 0) = 0.0372221;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(8333);
    tri3_xyze(0, 1) = 0.143815;
    tri3_xyze(1, 1) = -7.56737e-17;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = 0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-672);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143815;
    tri3_xyze(1, 0) = -7.56737e-17;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 1) = 0.140451;
    tri3_xyze(1, 1) = -9.23291e-17;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = 0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-672);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.140451;
    tri3_xyze(1, 0) = -9.23291e-17;
    tri3_xyze(2, 0) = 0.770611;
    nids.push_back(5581);
    tri3_xyze(0, 1) = 0.135665;
    tri3_xyze(1, 1) = 0.0363514;
    tri3_xyze(2, 1) = 0.770611;
    nids.push_back(8331);
    tri3_xyze(0, 2) = 0.139712;
    tri3_xyze(1, 2) = 0.0183934;
    tri3_xyze(2, 2) = 0.773262;
    nids.push_back(-672);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.146489;
    tri3_xyze(1, 0) = -5.78249e-17;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 1) = 0.143815;
    tri3_xyze(1, 1) = -7.56737e-17;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = 0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-673);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143815;
    tri3_xyze(1, 0) = -7.56737e-17;
    tri3_xyze(2, 0) = 0.775912;
    nids.push_back(5583);
    tri3_xyze(0, 1) = 0.138915;
    tri3_xyze(1, 1) = 0.0372221;
    tri3_xyze(2, 1) = 0.775912;
    nids.push_back(8333);
    tri3_xyze(0, 2) = 0.142679;
    tri3_xyze(1, 2) = 0.0187841;
    tri3_xyze(2, 2) = 0.778753;
    nids.push_back(-673);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.143372;
    tri3_xyze(1, 0) = 0.0384163;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(8337);
    tri3_xyze(0, 1) = 0.148429;
    tri3_xyze(1, 1) = -3.90641e-17;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = 0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-674);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.148429;
    tri3_xyze(1, 0) = -3.90641e-17;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 1) = 0.146489;
    tri3_xyze(1, 1) = -5.78249e-17;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = 0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-674);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.146489;
    tri3_xyze(1, 0) = -5.78249e-17;
    tri3_xyze(2, 0) = 0.781594;
    nids.push_back(5585);
    tri3_xyze(0, 1) = 0.141497;
    tri3_xyze(1, 1) = 0.0379141;
    tri3_xyze(2, 1) = 0.781594;
    nids.push_back(8335);
    tri3_xyze(0, 2) = 0.144947;
    tri3_xyze(1, 2) = 0.0190826;
    tri3_xyze(2, 2) = 0.78458;
    nids.push_back(-674);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144508;
    tri3_xyze(1, 0) = 0.0387208;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(8339);
    tri3_xyze(0, 1) = 0.149606;
    tri3_xyze(1, 1) = -1.96873e-17;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = 0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-675);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.149606;
    tri3_xyze(1, 0) = -1.96873e-17;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 1) = 0.148429;
    tri3_xyze(1, 1) = -3.90641e-17;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = 0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-675);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.148429;
    tri3_xyze(1, 0) = -3.90641e-17;
    tri3_xyze(2, 0) = 0.787566;
    nids.push_back(5587);
    tri3_xyze(0, 1) = 0.143372;
    tri3_xyze(1, 1) = 0.0384163;
    tri3_xyze(2, 1) = 0.787566;
    nids.push_back(8337);
    tri3_xyze(0, 2) = 0.146479;
    tri3_xyze(1, 2) = 0.0192843;
    tri3_xyze(2, 2) = 0.790649;
    nids.push_back(-675);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.144889;
    tri3_xyze(1, 0) = 0.0388229;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(8341);
    tri3_xyze(0, 1) = 0.15;
    tri3_xyze(1, 1) = 0;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = 0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-676);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.15;
    tri3_xyze(1, 0) = 0;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 1) = 0.149606;
    tri3_xyze(1, 1) = -1.96873e-17;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = 0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-676);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.149606;
    tri3_xyze(1, 0) = -1.96873e-17;
    tri3_xyze(2, 0) = 0.793733;
    nids.push_back(5589);
    tri3_xyze(0, 1) = 0.144508;
    tri3_xyze(1, 1) = 0.0387208;
    tri3_xyze(2, 1) = 0.793733;
    nids.push_back(8339);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = 0.0193859;
    tri3_xyze(2, 2) = 0.796867;
    nids.push_back(-676);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.149606;
    tri3_xyze(1, 0) = 1.96873e-17;
    tri3_xyze(2, 0) = 0.806267;
    nids.push_back(5593);
    tri3_xyze(0, 1) = 0.15;
    tri3_xyze(1, 1) = 0;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = 0.0193859;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-677);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = 0.15;
    tri3_xyze(1, 0) = 0;
    tri3_xyze(2, 0) = 0.8;
    nids.push_back(5591);
    tri3_xyze(0, 1) = 0.144889;
    tri3_xyze(1, 1) = 0.0388229;
    tri3_xyze(2, 1) = 0.8;
    nids.push_back(8341);
    tri3_xyze(0, 2) = 0.147251;
    tri3_xyze(1, 2) = 0.0193859;
    tri3_xyze(2, 2) = 0.803133;
    nids.push_back(-677);
    intersection.add_cut_side(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.2;
    hex8_xyze(1, 0) = -0.1;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1822);
    hex8_xyze(0, 1) = 0.2;
    hex8_xyze(1, 1) = -0.05;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1824);
    hex8_xyze(0, 2) = 0.15;
    hex8_xyze(1, 2) = -0.05;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1841);
    hex8_xyze(0, 3) = 0.15;
    hex8_xyze(1, 3) = -0.1;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1840);
    hex8_xyze(0, 4) = 0.2;
    hex8_xyze(1, 4) = -0.1;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1943);
    hex8_xyze(0, 5) = 0.2;
    hex8_xyze(1, 5) = -0.05;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1945);
    hex8_xyze(0, 6) = 0.15;
    hex8_xyze(1, 6) = -0.05;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1962);
    hex8_xyze(0, 7) = 0.15;
    hex8_xyze(1, 7) = -0.1;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1961);

    intersection.add_element(6889, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.2;
    hex8_xyze(1, 0) = -0.05;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1824);
    hex8_xyze(0, 1) = 0.2;
    hex8_xyze(1, 1) = 0;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1826);
    hex8_xyze(0, 2) = 0.15;
    hex8_xyze(1, 2) = -8.67362e-19;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1842);
    hex8_xyze(0, 3) = 0.15;
    hex8_xyze(1, 3) = -0.05;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1841);
    hex8_xyze(0, 4) = 0.2;
    hex8_xyze(1, 4) = -0.05;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1945);
    hex8_xyze(0, 5) = 0.2;
    hex8_xyze(1, 5) = 0;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1947);
    hex8_xyze(0, 6) = 0.15;
    hex8_xyze(1, 6) = -1.04083e-18;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1963);
    hex8_xyze(0, 7) = 0.15;
    hex8_xyze(1, 7) = -0.05;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1962);

    intersection.add_element(6890, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.2;
    hex8_xyze(1, 0) = 0;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1826);
    hex8_xyze(0, 1) = 0.2;
    hex8_xyze(1, 1) = 0.05;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1828);
    hex8_xyze(0, 2) = 0.15;
    hex8_xyze(1, 2) = 0.05;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1843);
    hex8_xyze(0, 3) = 0.15;
    hex8_xyze(1, 3) = -8.67362e-19;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1842);
    hex8_xyze(0, 4) = 0.2;
    hex8_xyze(1, 4) = 0;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1947);
    hex8_xyze(0, 5) = 0.2;
    hex8_xyze(1, 5) = 0.05;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1949);
    hex8_xyze(0, 6) = 0.15;
    hex8_xyze(1, 6) = 0.05;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1964);
    hex8_xyze(0, 7) = 0.15;
    hex8_xyze(1, 7) = -1.04083e-18;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1963);

    intersection.add_element(6891, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.15;
    hex8_xyze(1, 0) = -0.05;
    hex8_xyze(2, 0) = 0.75;
    nids.push_back(1841);
    hex8_xyze(0, 1) = 0.15;
    hex8_xyze(1, 1) = -8.67362e-19;
    hex8_xyze(2, 1) = 0.75;
    nids.push_back(1842);
    hex8_xyze(0, 2) = 0.1;
    hex8_xyze(1, 2) = 1.73472e-18;
    hex8_xyze(2, 2) = 0.75;
    nids.push_back(1853);
    hex8_xyze(0, 3) = 0.1;
    hex8_xyze(1, 3) = -0.05;
    hex8_xyze(2, 3) = 0.75;
    nids.push_back(1852);
    hex8_xyze(0, 4) = 0.15;
    hex8_xyze(1, 4) = -0.05;
    hex8_xyze(2, 4) = 0.8;
    nids.push_back(1962);
    hex8_xyze(0, 5) = 0.15;
    hex8_xyze(1, 5) = -1.04083e-18;
    hex8_xyze(2, 5) = 0.8;
    nids.push_back(1963);
    hex8_xyze(0, 6) = 0.1;
    hex8_xyze(1, 6) = 2.08167e-18;
    hex8_xyze(2, 6) = 0.8;
    nids.push_back(1974);
    hex8_xyze(0, 7) = 0.1;
    hex8_xyze(1, 7) = -0.05;
    hex8_xyze(2, 7) = 0.8;
    nids.push_back(1973);

    intersection.add_element(6900, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = 0.2;
    hex8_xyze(1, 0) = -0.05;
    hex8_xyze(2, 0) = 0.8;
    nids.push_back(1945);
    hex8_xyze(0, 1) = 0.2;
    hex8_xyze(1, 1) = 0;
    hex8_xyze(2, 1) = 0.8;
    nids.push_back(1947);
    hex8_xyze(0, 2) = 0.15;
    hex8_xyze(1, 2) = -1.04083e-18;
    hex8_xyze(2, 2) = 0.8;
    nids.push_back(1963);
    hex8_xyze(0, 3) = 0.15;
    hex8_xyze(1, 3) = -0.05;
    hex8_xyze(2, 3) = 0.8;
    nids.push_back(1962);
    hex8_xyze(0, 4) = 0.2;
    hex8_xyze(1, 4) = -0.05;
    hex8_xyze(2, 4) = 0.85;
    nids.push_back(2066);
    hex8_xyze(0, 5) = 0.2;
    hex8_xyze(1, 5) = 0;
    hex8_xyze(2, 5) = 0.85;
    nids.push_back(2068);
    hex8_xyze(0, 6) = 0.15;
    hex8_xyze(1, 6) = -1.21431e-18;
    hex8_xyze(2, 6) = 0.85;
    nids.push_back(2084);
    hex8_xyze(0, 7) = 0.15;
    hex8_xyze(1, 7) = -0.05;
    hex8_xyze(2, 7) = 0.85;
    nids.push_back(2083);

    intersection.add_element(6990, nids, hex8_xyze, Core::FE::CellType::hex8);
  }


  intersection.cut_test_cut(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation);
  intersection.cut_finalize(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation, false, true);
}
