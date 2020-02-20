TCut basic_cuts = "star_bbcwide > 0.08 && !gl1p_bad_run && !gl1p_bad_bunch "
                  "&& crossing != 1 && crossing < 111 ";

TCut run_qa_cuts_15em4 = "star_bbc30run / gl1_bbc30live > 1 && "
                         "star_bbc30run / gl1_bbc30live < 1.015 && "
                         "star_bbcwiderun / gl1_bbcwidelive > 1 && "
                         "star_bbcwiderun / gl1_bbcwidelive < 1.015 && "
                         "star_zdcwiderun / gl1_zdcwidelive > 1 && "
                         "star_zdcwiderun / gl1_zdcwidelive < 1.015  && "
                         "star_clkrun / gl1_clocklive > 1  && "
                         "star_clkrun / gl1_clocklive < 1.015 ";

TCut run_qa_cuts_10em4 = "star_bbc30run / gl1_bbc30live > 1 && "
                         "star_bbc30run / gl1_bbc30live < 1.010 && "
                         "star_bbcwiderun / gl1_bbcwidelive > 1 && "
                         "star_bbcwiderun / gl1_bbcwidelive < 1.010 && "
                         "star_zdcwiderun / gl1_zdcwidelive > 1 && "
                         "star_zdcwiderun / gl1_zdcwidelive < 1.010  && "
                         "star_clkrun / gl1_clocklive > 1  && "
                         "star_clkrun / gl1_clocklive < 1.010 ";

TCut run_qa_cuts_5em4 = "star_bbc30run / gl1_bbc30live > 1 && "
                        "star_bbc30run / gl1_bbc30live < 1.005 && "
                        "star_bbcwiderun / gl1_bbcwidelive > 1 && "
                        "star_bbcwiderun / gl1_bbcwidelive < 1.005 && "
                        "star_zdcwiderun / gl1_zdcwidelive > 1 && "
                        "star_zdcwiderun / gl1_zdcwidelive < 1.005  && "
                        "star_clkrun / gl1_clocklive > 1  && "
                        "star_clkrun / gl1_clocklive < 1.005 ";

TCut run_qa_cuts_1em4 = "star_bbc30run / gl1_bbc30live > 1 && "
                        "star_bbc30run / gl1_bbc30live < 1.001 && "
                        "star_bbcwiderun / gl1_bbcwidelive > 1 && "
                        "star_bbcwiderun / gl1_bbcwidelive < 1.001 && "
                        "star_zdcwiderun / gl1_zdcwidelive > 1 && "
                        "star_zdcwiderun / gl1_zdcwidelive < 1.001  && "
                        "star_clkrun / gl1_clocklive > 1  && "
                        "star_clkrun / gl1_clocklive < 1.001 ";

TCut run_qa_extra_runs =
    "runnumber != 386941 && runnumber != 387081 && runnumber != 387083 && "
    "runnumber != 387086 && runnumber != 387247 && runnumber != 387558 && "
    "runnumber != 389469 && runnumber != 390038 && runnumber != 390955 && "
    "runnumber != 391170 && runnumber != 391288 && runnumber != 391291 && "
    "runnumber != 391293 && runnumber != 391296 && runnumber != 391818 && "
    "runnumber != 392218 && runnumber != 392220 && runnumber != 392223 && "
    "runnumber != 392667 && runnumber != 392668 && runnumber != 393798 && "
    "runnumber != 394060 && runnumber != 394525 && runnumber != 394528 && "
    "runnumber != 394529 && runnumber != 394531 && runnumber != 395775 && "
    "runnumber != 396061 && runnumber != 396785 && runnumber != 397176 && "
    "runnumber != 397989 && runnumber != 397990";

// run_qa_extra_runs += "runnumber != 391175 && runnumber != 392281 &&"
//   "runnumber!=  392282 && runnumber != 392294 && runnumber 392296 &&
//   runnumber != 392299";

TCut crossing_qa_cuts = "star_zdc30cnt/gl1p_zdc_narrow > 0.9999 && "
                        "star_zdc30cnt/gl1p_zdc_narrow < 1.001 && "
                        "star_zdcwidecnt/gl1p_zdc_wide > 0.9999 && "
                        "star_zdcwidecnt/gl1p_zdc_wide < 1.001 && "
                        "star_bbc30cnt/gl1p_bbc_30 > 0.9999 && "
                        "star_bbc30cnt/gl1p_bbc_30 < 1.001 ";
// star_bbn:star_bbcwide
TCut single_arm_qa_cuts =
    "(star_bbn + 0.00596036 - 1.27234 * star_bbcwide + 0.27139 * star_bbcwide "
    "* star_bbcwide) < .008 &&"
    "(star_bbn + 0.00596036 - 1.27234 * star_bbcwide + 0.27139 * star_bbcwide "
    "* star_bbcwide) > -.008 &&"
    "(star_bbs < .0605 + 1.0625 * star_bbcwide) && star_zdn < 0.6";

TCut north_south_qa_cuts = "star_zdn / star_zds < .98 && star_zdn / star_zds > "
                           ".946 && star_bbn / star_bbs > .993 && star_bbn / "
                           "star_bbs < 1.005";

TCut default_cuts = basic_cuts + run_qa_cuts_15em4 + run_qa_extra_runs;
TCut tight_cuts = basic_cuts + run_qa_cuts_1em4 + run_qa_extra_runs;

TCut run_livetime_cuts =
    "gl1_clocklive/gl1_clockraw < 1 && gl1_clocklive/gl1_clockraw > 0.7 ";

TCut crossing_livetime_cuts = "star_clk/star_rclk > 0.7 && "
                              "star_bbcwidecnt/star_rbbcwidecnt > 0.7 && "
                              "star_zdcwidecnt/star_rzdcwidecnt > 0.7 && "
                              "star_bbc30cnt/star_rbbc30cnt > 0.7 && "
                              "star_zdc30cnt/star_rzdc30cnt > 0.7 && "
                              "star_bbscnt/star_rbbscnt > 0.7 && "
                              "star_bbcns2cnt/star_rbbcns2cnt > 0.7 && "
                              "star_zdscnt/star_rzdscnt > 0.7 && "
                              "star_zdncnt/star_rzdncnt > 0.7 && "
                              "star_zdcns2cnt/star_rzdcns2cnt > 0.7 ";

TCut bad_zdc_bbc_all = "runnumber != 397938 && runnumber != 397989 && "
                       "runnumber != 397000 && runnumber != 388838";

TCut pileup_cuts = "xbx_bbckn < .25 && xbx_bbckn > 0.215"
                   "&& xbx_bbcknerr > 0 && xbx_bbcknerr < 10"
                   "&& xbx_bbcks < .5 && xbx_bbcks > 0."
                   "&& xbx_zdckn< 10 && xbx_zdckn > 3.6"
                   "&& xbx_zdcks< 5 && xbx_zdcks > 1.5"
                   "&& xbx_zdckn > 3.78 - 0.23 * star_bbcwidecorr";

// TCut bad_fills = "gl1p_fill != 17399 && gl1p_fill != 17284 && gl1p_fill !=
// 17370 && gl1p_fill != 17415 && gl1p_fill != 17414"
//   "&& gl1p_fill != 17439 && gl1p_fill != 17407";
TCut bad_fills =
    "gl1p_fill != 17284 && gl1p_fill != 17293 && gl1p_fill != 17333 && "
    "gl1p_fill != 17338 && gl1p_fill != 17488"
    " && gl1p_fill != 17256 && gl1p_fill != 17261 && gl1p_fill != 17308 && "
    "gl1p_fill != 17359 && gl1p_fill != 17402 && gl1p_fill != 17403 && "
    "gl1p_fill != 17410 && gl1p_fill != 17415 && gl1p_fill != 17417 && "
    "gl1p_fill != 17545 && gl1p_fill != 17534 && gl1p_fill != 17370 && "
    "gl1p_fill != 17399 && gl1p_fill != 17519 && gl1p_fill != 17346 && "
    "gl1p_fill != 17341";

TCut fill_cuts = "!(gl1p_fill == 17389 && star_bbcwidecorr < .3)"
                 " && !(gl1p_fill == 17423 && crossing == 13)"
                 " && !(gl1p_fill == 17423 && crossing >= 31 && crossing <= 39)"
                 " && !(gl1p_fill == 17430 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17452 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17453 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17454 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17461 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17466 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17466 && star_bbcwidecorr < .33)"
                 " && !(gl1p_fill == 17467 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17467 && star_bbcwidecorr > .3)"
                 " && !(gl1p_fill == 17467 && xbx_bbckn > .235)"
                 " && !(gl1p_fill == 17472 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17473 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17482 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17484 && crossing >= 71 && crossing <= 79)"
                 " && !(gl1p_fill == 17512 && star_bbcwidecorr < .43)"
                 " && !(gl1p_fill == 17514 && star_bbcwidecorr < .3)"
                 " && !(gl1p_fill == 17515 && star_bbcwidecorr < .4)"
                 " && !(gl1p_fill == 17517 && star_bbcwidecorr < .3)"
                 " && !(gl1p_fill == 17518 && star_bbcwidecorr < .3)"
                 " && !(gl1p_fill == 17519 && star_bbcwidecorr < .3)"
                 " && fbf_bbckn != -999";

TCut one_run_fills = "gl1p_fill != 17222"
                     "&& gl1p_fill != 17223"
                     "&& gl1p_fill != 17261"
                     "&& gl1p_fill != 17328"
                     "&& gl1p_fill != 17394"
                     "&& gl1p_fill != 17403"
                     "&& gl1p_fill != 17438"
                     "&& gl1p_fill != 17461"
                     "&& gl1p_fill != 17482"
                     "&& gl1p_fill != 17503"
                     "&& gl1p_fill != 17530"
                     "&& gl1p_fill != 17536";

TCut two_run_fills = one_run_fills + "gl1p_fill != 17256"
                                     "&& gl1p_fill != 17281"
                                     "&& gl1p_fill != 17346"
                                     "&& gl1p_fill != 17389"
                                     "&& gl1p_fill != 17439"
                                     "&& gl1p_fill != 17467"
                                     "&& gl1p_fill != 17492"
                                     "&& gl1p_fill != 17515"
                                     "&& gl1p_fill != 17518"
                                     "&& gl1p_fill != 17570"
                                     "&& gl1p_fill != 17571";

TCut three_run_fills = two_run_fills + "gl1p_fill != 17219"
                                       "&& gl1p_fill != 17253"
                                       "&& gl1p_fill != 17297"
                                       "&& gl1p_fill != 17317"
                                       "&& gl1p_fill != 17396"
                                       "&& gl1p_fill != 17407"
                                       "&& gl1p_fill != 17451"
                                       "&& gl1p_fill != 17452"
                                       "&& gl1p_fill != 17453"
                                       "&& gl1p_fill != 17466"
                                       "&& gl1p_fill != 17470"
                                       "&& gl1p_fill != 17479"
                                       "&& gl1p_fill != 17502"
                                       "&& gl1p_fill != 17513"
                                       "&& gl1p_fill != 17514"
                                       "&& gl1p_fill != 17524"
                                       "&& gl1p_fill != 17527"
                                       "&& gl1p_fill != 17529"
                                       "&& gl1p_fill != 17560"
                                       "&& gl1p_fill != 17593";

TCut four_run_fills = three_run_fills + "gl1p_fill != 17236"
                                        "&& gl1p_fill != 17240"
                                        "&& gl1p_fill != 17247"
                                        "&& gl1p_fill != 17306"
                                        "&& gl1p_fill != 17308"
                                        "&& gl1p_fill != 17311"
                                        "&& gl1p_fill != 17322"
                                        "&& gl1p_fill != 17331"
                                        "&& gl1p_fill != 17380"
                                        "&& gl1p_fill != 17410"
                                        "&& gl1p_fill != 17423"
                                        "&& gl1p_fill != 17427"
                                        "&& gl1p_fill != 17434"
                                        "&& gl1p_fill != 17436"
                                        "&& gl1p_fill != 17484"
                                        "&& gl1p_fill != 17488"
                                        "&& gl1p_fill != 17491"
                                        "&& gl1p_fill != 17512"
                                        "&& gl1p_fill != 17538"
                                        "&& gl1p_fill != 17544"
                                        "&& gl1p_fill != 17548"
                                        "&& gl1p_fill != 17558"
                                        "&& gl1p_fill != 17587"
                                        "&& gl1p_fill != 17594";

TCut before_spinpats = "gl1p_fill >=  17256";

TCut bad_bunch_shuffling =
    "bsh_factor < 1.2"; //"runnumber != 388548 && runnumber != 389904 &&
                        //runnumber != 390026 && runnumber != 390029 &&
                        //runnumber != 390030 && runnumber != 390032 &&
                        //runnumber != 390039 && runnumber != 395402";

TCut final_cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
                  crossing_livetime_cuts + single_arm_qa_cuts +
                  north_south_qa_cuts + pileup_cuts + bad_fills + fill_cuts +
                  before_spinpats + bad_bunch_shuffling;

const int N_RUNS = 1087;
const int N_CROSSINGS = 120;

double ppmm[] = {+1, +1, -1, -1, +1, +1, -1, -1};
double mmpp[] = {-1, -1, +1, +1, -1, -1, +1, +1};
double pppp[] = {+1, +1, +1, +1, -1, -1, -1, -1};
double mmmm[] = {-1, -1, -1, -1, +1, +1, +1, +1};
double pppps2[] = {+1, +1, -1, -1, -1, -1, +1, +1};
double mmmms2[] = {-1, -1, +1, +1, +1, +1, -1, -1};

bool equal(double *v, double *a) {
  for (std::size_t i = 0; i < 8; i++) {
    // std::cout << "*(v + i) " << *(v + i) << std::endl;
    // std::cout << "*(a + i) " << *(a + i) << std::endl;
    if (*(v + i) != *(a + i))
      return false;
  }
  return true;
}

struct rl_db {
  TFile *f;
  TTree *t;
  //there is also an int m[] defined in protected, below.
  rl_db() {
    f = TFile::Open("/phenix/spin2/pmontu/offline/analysis/pmontu/"
                    "relative_luminosity/SpinDB/unique_db.root",
                    "read");
    t = (TTree *)f->Get("t");
    load_map();
  }
  int get_tree_entry(const int run_index, const int crossing) {
    return run_index * 120 + crossing;
  }

  int get_spin_pattern(const int runnumber) {
    int result = 0;
    int index = m[runnumber];
    int entry = index * 120;
    // std::cout << "index " << index << std::endl;
    // std::cout << "entry " << entry << std::endl;

    //load the 8 events starting at event 'entry' into t's temporary vectors:
    int n_points = t->Draw("bpat:ypat", "", "", 8, entry);
    //now V1 contains eight values of bpat and V2 contains 8 values of ypat.

    // std::cout << "*(t->GetV1() + 3) " << *(t->GetV1() + 3) << std::endl;
    int ishift = 2;
    if (equal(t->GetV1(), ppmm) && equal(t->GetV2(), pppp)) {
      // cout << runnumber << " pat1 " << endl;
      return 1;
    } else if (equal(t->GetV1(), mmpp) && equal(t->GetV2(), pppp)) {
      // cout << runnumber << " pat2 " << endl;
      return 2;
    } else if (equal(t->GetV1(), ppmm) && equal(t->GetV2(), mmmm)) {
      // cout << runnumber << " pat3 " << endl;
      return 3;
    } else if (equal(t->GetV1(), mmpp) && equal(t->GetV2(), mmmm)) {
      // cout << runnumber << " pat4 " << endl;
      return 4;
    } else if (equal(t->GetV1(), pppp) && equal(t->GetV2(), ppmm)) {
      // cout << runnumber << " pat5 " << endl;
      return 5;
    } else if (equal(t->GetV1(), pppp) && equal(t->GetV2(), mmpp)) {
      // cout << runnumber << " pat6 " << endl;
      return 6;
    } else if (equal(t->GetV1(), mmmm) && equal(t->GetV2(), ppmm)) {
      // cout << runnumber << " pat7 " << endl;
      return 7;
    } else if (equal(t->GetV1(), mmmm) && equal(t->GetV2(), mmpp)) {
      // cout << runnumber << " pat8 " << endl;
      return 8;
    } else if (equal(t->GetV1(), mmpp) && equal(t->GetV2(), pppps2)) {
      // cout << runnumber << " pat21 " << endl;
      return 21;
    } else if (equal(t->GetV1(), ppmm) && equal(t->GetV2(), pppps2)) {
      // cout << runnumber << " pat22 " << endl;
      return 22;
    } else if (equal(t->GetV1(), mmpp) && equal(t->GetV2(), mmmms2)) {
      // cout << runnumber << " pat23 " << endl;
      return 23;
    } else if (equal(t->GetV1(), ppmm) && equal(t->GetV2(), mmmms2)) {
      // cout << runnumber << " pat24 " << endl;
      return 24;
    } else if (equal(t->GetV1(), pppps2) && equal(t->GetV2(), mmpp)) {
      // cout << runnumber << " pat25 " << endl;
      return 25;
    } else if (equal(t->GetV1(), pppps2) && equal(t->GetV2(), ppmm)) {
      // cout << runnumber << " pat26 " << endl;
      return 26;
    } else if (equal(t->GetV1(), mmmms2) && equal(t->GetV2(), mmpp)) {
      // cout << runnumber << " pat27 " << endl;
      return 27;
    } else if (equal(t->GetV1(), mmmms2) && equal(t->GetV2(), ppmm)) {
      // cout << runnumber << " pat28 " << endl;
      return 28;
    } else {
      // cout << "NO MATCH!!\n";
      return -999;
    }
  }

  TGraphErrors *get_ratio_vs_crossing(const int runnumber, const TString num,
                                      const TString den,
                                      const TCut extra_cut = "",
                                      const double err_factor = 1.) {

    // TString run_str = "runnumber == ";
    // run_str += runnumber;
    // TCut run_cut = run_str;

    TCut cuts = final_cuts + extra_cut;

    int index = m[runnumber];

    TString nume = num;
    nume += "err * ";
    nume += err_factor;
    TString dene = den;
    dene += "err * ";
    dene += err_factor;

    TString formula = num + " / " + den + ":" + nume + "/" + num + ":" + dene +
                      "/" + den + ":crossing";
    // std::cout << "formula " << formula << std::endl;
    int n_p = t->Draw(formula, cuts, "", 120, index * 120);
    // std::cout << "n_p " << n_p << std::endl;

    if (n_p == 0)
      return new TGraphErrors(0);

    vector<double> err;

    for (std::size_t ip = 0; ip < n_p; ip++) {
      double errnum = t->GetV2()[ip];
      double errden = t->GetV3()[ip];
      err.push_back(fabs(t->GetV1()[ip]) *
                    sqrt(pow(errnum, 2) + pow(errden, 2)));
      // std::cout << "ratio " << t->GetV1()[ip] << std::endl;
      // std::cout << "errnum " << errnum << std::endl;
      // std::cout << "errden " << errden << std::endl;
      // std::cout << "err[ip] " << err[ip] << std::endl;
    }

    vector<double> zero_pad(120, 0);

    return new TGraphErrors(n_p, t->GetV4(), t->GetV1(), &zero_pad[0], &err[0]);
  }

  TGraphErrors *get_al(const bool is_blue = true,
                       const TString num = "fbf_bbcwide",
                       const TString den = "fbf_zdcwide") {
    TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
                crossing_livetime_cuts + single_arm_qa_cuts +
                north_south_qa_cuts + pileup_cuts + bad_fills + fill_cuts +
                before_spinpats + bad_bunch_shuffling;

    TString pat = (is_blue ? "bpat" : "ypat");

    TCut pcuts = cuts && TString(pat + "==1");
    TCut mcuts = cuts && TString(pat + "==-1");

    int npp = t->Draw("runnumber", final_cuts + "bpat*ypat==1");
    vector<int> runspp;
    for (std::size_t i = 0; i < npp; i++) {
      runspp.push_back(t->GetV1()[i]);
    }

    int npm = t->Draw("runnumber", final_cuts + "bpat*ypat==-1");
    vector<int> runspm;
    for (std::size_t i = 0; i < npm; i++) {
      runspm.push_back(t->GetV1()[i]);
    }

    TGraphErrors *result = new TGraphErrors();

    TString bbc_scaler = num;
    bbc_scaler += " * star_clk";
    TString zdc_scaler = den;
    zdc_scaler += " * star_clk";

    TString bbc_scaler_err = num;
    bbc_scaler_err += "err * star_clk";
    TString zdc_scaler_err = den;
    zdc_scaler_err += "err * star_clk";

    // std::cout << "bbc_scaler " << bbc_scaler << std::endl;
    // std::cout << "zdc_scaler " << zdc_scaler << std::endl;
    // std::cout << "bbc_scaler_err " << bbc_scaler_err << std::endl;
    // std::cout << "zdc_scaler_err " << zdc_scaler_err << std::endl;

    vector<double> runs_p;
    vector<double> cros_p;
    vector<double> zdc_p;
    vector<double> zdc_p_err;
    vector<double> bbc_p;
    vector<double> bbc_p_err;
    vector<double> runs_m;
    vector<double> cros_m;
    vector<double> zdc_m;
    vector<double> zdc_m_err;
    vector<double> bbc_m;
    vector<double> bbc_m_err;

    TString query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

    int n_points = t->Draw(query, pcuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      runs_p.push_back(t->GetV1()[ip]);
      cros_p.push_back(t->GetV2()[ip]);
      zdc_p.push_back(t->GetV3()[ip]);
      zdc_p_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    int n_points2 = t->Draw(query, pcuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      if (t->GetV1()[ip] != runs_p[ip])
        cout << "AAA\n";
      if (t->GetV2()[ip] != cros_p[ip])
        cout << "AAA\n";
      bbc_p.push_back(t->GetV3()[ip]);
      bbc_p_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

    n_points = t->Draw(query, mcuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      runs_m.push_back(t->GetV1()[ip]);
      cros_m.push_back(t->GetV2()[ip]);
      zdc_m.push_back(t->GetV3()[ip]);
      zdc_m_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    n_points2 = t->Draw(query, mcuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      if (t->GetV1()[ip] != runs_m[ip])
        cout << "AAA\n";
      if (t->GetV2()[ip] != cros_m[ip])
        cout << "AAA\n";
      bbc_m.push_back(t->GetV3()[ip]);
      bbc_m_err.push_back(t->GetV4()[ip]);
    }

    // std::size_t index_pp = 0;
    // std::size_t index_m = 0;

    map<int, int>::iterator rit = m.begin();

    while (rit != m.end()) {

      int cnpp = 0;
      int cnpm = 0;

      int currun = rit->first;

      for (std::size_t i = 0; i < npp; i++) {
        if (runspp[i] == currun) {
          cnpp++;
        }
      }

      for (std::size_t i = 0; i < npm; i++) {
        if (runspm[i] == currun) {
          cnpm++;
        }
      }
      if (cnpp > 45 && cnpm > 45) {
        // if (currun > 386946) break;
        TGraphErrors zplus;
        TGraphErrors zminu;
        TGraphErrors bplus;
        TGraphErrors bminu;

        for (std::size_t ip = 0; ip < runs_p.size(); ip++) {
          if (runs_p[ip] == currun) {
            zplus.SetPoint(zplus.GetN(), cros_p[ip], zdc_p[ip]);
            zplus.SetPointError(zplus.GetN() - 1, 0, zdc_p_err[ip]);
            bplus.SetPoint(bplus.GetN(), cros_p[ip], bbc_p[ip]);
            bplus.SetPointError(bplus.GetN() - 1, 0, bbc_p_err[ip]);
          }
        }

        for (std::size_t ip = 0; ip < runs_m.size(); ip++) {
          if (runs_m[ip] == currun) {
            zminu.SetPoint(zminu.GetN(), cros_m[ip], zdc_m[ip]);
            zminu.SetPointError(zminu.GetN() - 1, 0, zdc_m_err[ip]);
            bminu.SetPoint(bminu.GetN(), cros_m[ip], bbc_m[ip]);
            bminu.SetPointError(bminu.GetN() - 1, 0, bbc_m_err[ip]);
          }
        }

        if (zplus.GetN() > 41 && zminu.GetN() > 43) {

          zplus.Fit("pol0", "q");
          zminu.Fit("pol0", "q");
          bplus.Fit("pol0", "q");
          bminu.Fit("pol0", "q");

          double dzdc_p = 0.;
          double dzdc_p_err_sq = 0.;
          double dzdc_m = 0.;
          double dzdc_m_err_sq = 0.;
          double dbbc_p = 0.;
          double dbbc_p_err_sq = 0.;
          double dbbc_m = 0.;
          double dbbc_m_err_sq = 0.;

          for (std::size_t ip = 0; ip < zplus.GetN(); ip++) {
            dzdc_p += zplus.GetY()[ip];
            dzdc_p_err_sq += pow(zplus.GetEY()[ip], 2);
            dbbc_p += bplus.GetY()[ip];
            dbbc_p_err_sq += pow(bplus.GetEY()[ip], 2);
          }

          for (std::size_t ip = 0; ip < zminu.GetN(); ip++) {
            dzdc_m += zminu.GetY()[ip];
            dzdc_m_err_sq += pow(zminu.GetEY()[ip], 2);
            dbbc_m += bminu.GetY()[ip];
            dbbc_m_err_sq += pow(bminu.GetEY()[ip], 2);
          }

          double dzdc_p_err = sqrt(dzdc_p_err_sq);
          double dzdc_m_err = sqrt(dzdc_m_err_sq);
          double dbbc_p_err = sqrt(dbbc_p_err_sq);
          double dbbc_m_err = sqrt(dbbc_m_err_sq);

          TString beam = is_blue ? "gl1p_bpol" : "gl1p_ypol";
          TString beamerr = is_blue ? "gl1p_bpolerr" : "gl1p_ypolerr";

          t->GetEntry(rit->second * 120);

          double pol = t->GetLeaf(beam)->GetValue(0);
          double pol_err = t->GetLeaf(beamerr)->GetValue(0);

          double a_l = ((dzdc_p / dbbc_p) - (dzdc_m / dbbc_m)) /
                       (pol * ((dzdc_p / dbbc_p) + (dzdc_m / dbbc_m)));

          double pol_rel_error_sq = pow(pol_err / pol, 2);
          double scaler_rel_error_sq =
              pow(dzdc_p_err / dzdc_p, 2) + pow(dzdc_m_err / dzdc_m, 2) +
              pow(dbbc_p_err / dbbc_p, 2) + pow(dbbc_m_err / dbbc_m, 2);

          double radicant = pol_rel_error_sq;
          radicant +=
              pow(2 * dzdc_p * dzdc_m * dbbc_p * dbbc_m, 2) /
              pow(pow(dzdc_p * dbbc_m, 2) - pow(dzdc_m * dbbc_p, 2), 2) *
              scaler_rel_error_sq;
          double a_l_err = fabs(a_l) * sqrt(radicant);

          if (false) {
            std::cout << "dzdc_p " << zdc_p << std::endl;
            std::cout << "dzdc_p_err " << zdc_p_err << std::endl;
            std::cout << "dzdc_m " << zdc_m << std::endl;
            std::cout << "dzdc_m_err " << zdc_m_err << std::endl;
            std::cout << "dbbc_p " << bbc_p << std::endl;
            std::cout << "dbbc_p_err " << bbc_p_err << std::endl;
            std::cout << "dbbc_m " << bbc_m << std::endl;
            std::cout << "dbbc_m_err " << bbc_m_err << std::endl;
            std::cout << "pol_blue " << pol_blue << std::endl;
            std::cout << "pol_yell " << pol_yell << std::endl;
            std::cout << "pol_blue_err " << pol_blue_err << std::endl;
            std::cout << "pol_yell_err " << pol_yell_err << std::endl;
            std::cout << "a_ll " << a_ll << std::endl;
            std::cout << "pol_rel_error_sq " << pol_rel_error_sq << std::endl;
            std::cout << "scaler_rel_error_sq " << scaler_rel_error_sq
                      << std::endl;
            std::cout << "radicant " << radicant << std::endl;
            std::cout << "a_ll_err " << a_ll_err << std::endl;
          }
          if (true) {
            result->SetPoint(result->GetN(), currun, a_l);
            result->SetPointError(result->GetN() - 1, 0, a_l_err);
          }
        }
      }
      ++rit;
    }
    return result;
  }

  TGraphErrors *get_all(TCut extra_cut = "") {
    TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
                crossing_livetime_cuts + single_arm_qa_cuts +
                north_south_qa_cuts + extra_cut + pileup_cuts + bad_fills +
                fill_cuts + before_spinpats + bad_bunch_shuffling;

    TCut ppcuts = cuts + "bpat*ypat==1";
    TCut pmcuts = cuts + "bpat*ypat==-1";

    int npp = t->Draw("runnumber", final_cuts + "bpat*ypat==1");
    vector<int> runspp;
    for (std::size_t i = 0; i < npp; i++) {
      runspp.push_back(t->GetV1()[i]);
    }

    int npm = t->Draw("runnumber", final_cuts + "bpat*ypat==-1");
    vector<int> runspm;
    for (std::size_t i = 0; i < npm; i++) {
      runspm.push_back(t->GetV1()[i]);
    }

    TGraphErrors *result = new TGraphErrors();

    // TString bbc_scaler = "star_bbcwidecorr * star_clk";
    // TString zdc_scaler = "star_zdcwidecorr * star_clk";
    // TString bbc_scaler_err = "star_bbcwidecorrerr * star_clk";
    // TString zdc_scaler_err = "star_zdcwidecorrerr * star_clk";

    TString bbc_scaler = "fbf_bbcwide * star_clk";
    TString zdc_scaler = "fbf_zdcwide * star_clk";
    // TString bbc_scaler_err = "glob_bbcwideerr * star_clk";
    // TString zdc_scaler_err = "glob_zdcwideerr * star_clk";
    TString bbc_scaler_err = "fbf_bbcwideerr * star_clk";
    TString zdc_scaler_err = "fbf_zdcwideerr * star_clk";
    // TString bbc_scaler_err = "sqrt(fbf_bbcwideerr*fbf_bbcwideerr +
    // xbx_bbcwideerr*xbx_bbcwideerr) * star_clk";
    // TString zdc_scaler_err = "sqrt(fbf_zdcwideerr*fbf_zdcwideerr +
    // xbx_zdcwideerr*xbx_zdcwideerr) * star_clk";
    // TString bbc_scaler_err = "star_bbcwidecorrerr * star_clk";
    // TString zdc_scaler_err = "star_zdcwidecorrerr * star_clk";

    // TString bbc_scaler = "rbr_bbcwide * star_clk";
    // TString zdc_scaler = "rbr_zdcwide * star_clk";
    // TString bbc_scaler_err = "star_bbcwidecorrerr * star_clk";
    // TString zdc_scaler_err = "star_zdcwidecorrerr * star_clk";

    // TString bbc_scaler = "glob_bbcwide * star_clk";
    // TString zdc_scaler = "glob_zdcwide * star_clk";
    // TString bbc_scaler_err = "star_bbcwidecorrerr * star_clk";
    // TString zdc_scaler_err = "star_zdcwidecorrerr * star_clk";

    // TString bbc_scaler = "star_bbcwide * star_clk";
    // TString zdc_scaler = "star_zdcwide * star_clk";
    // TString bbc_scaler_err = "star_bbcwideerr * star_clk";
    // TString zdc_scaler_err = "star_zdcwideerr * star_clk";

    // TString bbc_scaler = "xbx_bbcwide * star_clk";
    // TString zdc_scaler = "xbx_zdcwide * star_clk";
    // TString bbc_scaler_err = "xbx_bbcwideerr * star_clk";
    // TString zdc_scaler_err = "xbx_zdcwideerr * star_clk";

    vector<double> runs_pp;
    vector<double> cros_pp;
    vector<double> zdc_pp;
    vector<double> zdc_pp_err;
    vector<double> bbc_pp;
    vector<double> bbc_pp_err;
    vector<double> runs_pm;
    vector<double> cros_pm;
    vector<double> zdc_pm;
    vector<double> zdc_pm_err;
    vector<double> bbc_pm;
    vector<double> bbc_pm_err;

    TString query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

    int n_points = t->Draw(query, ppcuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      runs_pp.push_back(t->GetV1()[ip]);
      cros_pp.push_back(t->GetV2()[ip]);
      zdc_pp.push_back(t->GetV3()[ip]);
      zdc_pp_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    int n_points2 = t->Draw(query, ppcuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      if (t->GetV1()[ip] != runs_pp[ip])
        cout << "AAA\n";
      if (t->GetV2()[ip] != cros_pp[ip])
        cout << "AAA\n";
      bbc_pp.push_back(t->GetV3()[ip]);
      bbc_pp_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

    n_points = t->Draw(query, pmcuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      runs_pm.push_back(t->GetV1()[ip]);
      cros_pm.push_back(t->GetV2()[ip]);
      zdc_pm.push_back(t->GetV3()[ip]);
      zdc_pm_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    n_points2 = t->Draw(query, pmcuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      if (t->GetV1()[ip] != runs_pm[ip])
        cout << "AAA\n";
      if (t->GetV2()[ip] != cros_pm[ip])
        cout << "AAA\n";
      bbc_pm.push_back(t->GetV3()[ip]);
      bbc_pm_err.push_back(t->GetV4()[ip]);
    }

    // std::size_t index_pp = 0;
    // std::size_t index_pm = 0;

    map<int, int>::iterator rit = m.begin();

    while (rit != m.end()) {
      int currun = rit->first;
      // if (currun > 386946) break;
      int cnpp = 0;
      int cnpm = 0;

      TGraphErrors zsame;
      TGraphErrors zoppo;
      TGraphErrors bsame;
      TGraphErrors boppo;

      for (std::size_t ip = 0; ip < runs_pp.size(); ip++) {
        if (runs_pp[ip] == currun) {
          zsame.SetPoint(zsame.GetN(), cros_pp[ip], zdc_pp[ip]);
          zsame.SetPointError(zsame.GetN() - 1, 0, zdc_pp_err[ip]);
          bsame.SetPoint(bsame.GetN(), cros_pp[ip], bbc_pp[ip]);
          bsame.SetPointError(bsame.GetN() - 1, 0, bbc_pp_err[ip]);
        }
      }

      for (std::size_t ip = 0; ip < runs_pm.size(); ip++) {
        if (runs_pm[ip] == currun) {
          zoppo.SetPoint(zoppo.GetN(), cros_pm[ip], zdc_pm[ip]);
          zoppo.SetPointError(zoppo.GetN() - 1, 0, zdc_pm_err[ip]);
          boppo.SetPoint(boppo.GetN(), cros_pm[ip], bbc_pm[ip]);
          boppo.SetPointError(boppo.GetN() - 1, 0, bbc_pm_err[ip]);
        }
      }

      for (std::size_t i = 0; i < npp; i++) {
        if (runspp[i] == currun) {
          cnpp++;
        }
      }

      for (std::size_t i = 0; i < npm; i++) {
        if (runspm[i] == currun) {
          cnpm++;
        }
      }

      // std::cout << "currun " << currun << " ";
      // cout << zsame.GetN() << " " << zoppo.GetN() << " " << bsame.GetN() << "
      // " << boppo.GetN() << endl;

      // std::cout << "cnpp " << cnpp << std::endl;
      // std::cout << "cnpm " << cnpm << std::endl;

      if (cnpp > 45 && cnpm > 45) {
        // if (zsame.GetN() > 45 && zoppo.GetN() > 45) {

        // TString cname =
        // "/phenix/WWW/p/draft/pmontu/relative_luminosity/hel_plots/";
        // cname += currun;
        // cname += ".png";
        // TCanvas *c2 = new TCanvas("c2", "c2", 800 * 2, 600 * 2);
        // c2->Divide(2, 2);
        zsame.Fit("pol0", "q");
        zoppo.Fit("pol0", "q");
        bsame.Fit("pol0", "q");
        boppo.Fit("pol0", "q");

        // c2->cd(1);
        // zsame.Draw("ape");
        // c2->cd(2);
        // zoppo.Draw("ape");
        // c2->cd(3);
        // bsame.Draw("ape");
        // c2->cd(4);
        // boppo.Draw("ape");

        // c2->SaveAs(cname);
        // delete c2;

        double zdc_sa = 0.;
        double zdc_sa_err_sq = 0.;
        double zdc_op = 0.;
        double zdc_op_err_sq = 0.;
        double bbc_sa = 0.;
        double bbc_sa_err_sq = 0.;
        double bbc_op = 0.;
        double bbc_op_err_sq = 0.;

        for (std::size_t ip = 0; ip < zsame.GetN(); ip++) {
          zdc_sa += zsame.GetY()[ip];
          zdc_sa_err_sq += pow(zsame.GetEY()[ip], 2);
          bbc_sa += bsame.GetY()[ip];
          bbc_sa_err_sq += pow(bsame.GetEY()[ip], 2);
        }

        for (std::size_t ip = 0; ip < zoppo.GetN(); ip++) {
          zdc_op += zoppo.GetY()[ip];
          zdc_op_err_sq += pow(zoppo.GetEY()[ip], 2);
          bbc_op += boppo.GetY()[ip];
          bbc_op_err_sq += pow(boppo.GetEY()[ip], 2);
        }

        double zdc_sa_err = sqrt(zdc_sa_err_sq);
        double zdc_op_err = sqrt(zdc_op_err_sq);
        double bbc_sa_err = sqrt(bbc_sa_err_sq);
        double bbc_op_err = sqrt(bbc_op_err_sq);

        // TF1 *fzsame = zsame.GetFunction("pol0");
        // TF1 *fzoppo = zoppo.GetFunction("pol0");
        // TF1 *fbsame = bsame.GetFunction("pol0");
        // TF1 *fboppo = boppo.GetFunction("pol0");

        // double zdc_sa =     fzsame->GetParameter(0);
        // double zdc_sa_err = fzsame->GetParError(0);
        // double zdc_op =     fzoppo->GetParameter(0);
        // double zdc_op_err = fzoppo->GetParError(0);
        // double bbc_sa =     fbsame->GetParameter(0);
        // double bbc_sa_err = fbsame->GetParError(0);
        // double bbc_op =     fboppo->GetParameter(0);
        // double bbc_op_err = fboppo->GetParError(0);

        t->GetEntry(rit->second * 120);
        double pol_blue = t->GetLeaf("gl1p_bpol")->GetValue(0);
        double pol_yell = t->GetLeaf("gl1p_ypol")->GetValue(0);
        double pol_blue_err = t->GetLeaf("gl1p_bpolerr")->GetValue(0);
        double pol_yell_err = t->GetLeaf("gl1p_ypolerr")->GetValue(0);

        double a_ll =
            ((zdc_sa / bbc_sa) - (zdc_op / bbc_op)) /
            (pol_blue * pol_yell * ((zdc_sa / bbc_sa) + (zdc_op / bbc_op)));

        double pol_rel_error_sq =
            pow(pol_blue_err / pol_blue, 2) + pow(pol_yell_err / pol_yell, 2);
        double scaler_rel_error_sq =
            pow(zdc_sa_err / zdc_sa, 2) + pow(zdc_op_err / zdc_op, 2) +
            pow(bbc_sa_err / bbc_sa, 2) + pow(bbc_op_err / bbc_op, 2);

        double radicant = pol_rel_error_sq;
        radicant += pow(2 * zdc_sa * zdc_op * bbc_sa * bbc_op, 2) /
                    pow(pow(zdc_sa * bbc_op, 2) - pow(zdc_op * bbc_sa, 2), 2) *
                    scaler_rel_error_sq;
        double a_ll_err = fabs(a_ll) * sqrt(radicant);

        if (false) {
          std::cout << "zdc_sa " << zdc_sa << std::endl;
          std::cout << "zdc_sa_err " << zdc_sa_err << std::endl;
          std::cout << "zdc_op " << zdc_op << std::endl;
          std::cout << "zdc_op_err " << zdc_op_err << std::endl;
          std::cout << "bbc_sa " << bbc_sa << std::endl;
          std::cout << "bbc_sa_err " << bbc_sa_err << std::endl;
          std::cout << "bbc_op " << bbc_op << std::endl;
          std::cout << "bbc_op_err " << bbc_op_err << std::endl;
          std::cout << "pol_blue " << pol_blue << std::endl;
          std::cout << "pol_yell " << pol_yell << std::endl;
          std::cout << "pol_blue_err " << pol_blue_err << std::endl;
          std::cout << "pol_yell_err " << pol_yell_err << std::endl;
          std::cout << "a_ll " << a_ll << std::endl;
          std::cout << "pol_rel_error_sq " << pol_rel_error_sq << std::endl;
          std::cout << "scaler_rel_error_sq " << scaler_rel_error_sq
                    << std::endl;
          std::cout << "radicant " << radicant << std::endl;
          std::cout << "a_ll_err " << a_ll_err << std::endl;
        }
        if (a_ll > -1 && a_ll < 1 && a_ll_err > 0 && a_ll_err < 1) {
          result->SetPoint(result->GetN(), currun, a_ll);
          result->SetPointError(result->GetN() - 1, 0, a_ll_err);
        }
      }

      ++rit;
    }
    return result;
  }

  TGraphErrors *get_syst_err_all(TCut extra_cut = "") {

    TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
                crossing_livetime_cuts + single_arm_qa_cuts +
                north_south_qa_cuts + extra_cut + pileup_cuts + bad_fills +
                fill_cuts + before_spinpats;

    TGraphErrors *result = new TGraphErrors();

    int itn = 0;

    for (std::size_t ib = 0; ib < 9; ib += 1) {
      for (std::size_t iz = 0; iz < 9; iz += 1) {

        itn++;

        cout << itn << "/81" << endl;

        TGraphErrors *currresult = new TGraphErrors();
        TString gname = "g";
        gname += itn;
        currresult->SetName(gname);

        TCut ppcuts = cuts + "bpat*ypat==1";
        TCut pmcuts = cuts + "bpat*ypat==-1";

        TString bbc_scaler = "fbf_bbcwide";
        bbc_scaler += ib;
        bbc_scaler += " * star_clk";
        TString zdc_scaler = "fbf_zdcwide";
        zdc_scaler += iz;
        zdc_scaler += " * star_clk";

        TString bbc_scaler_err = "fbf_bbcwideerr * star_clk";
        TString zdc_scaler_err = "fbf_zdcwideerr * star_clk";

        // std::cout << "bbc_scaler " << bbc_scaler << std::endl;
        // std::cout << "zdc_scaler " << zdc_scaler << std::endl;
        // std::cout << "curresult->GetName() " << curresult->GetName() <<
        // std::endl;

        vector<double> runs_pp;
        vector<double> cros_pp;
        vector<double> zdc_pp;
        vector<double> zdc_pp_err;
        vector<double> bbc_pp;
        vector<double> bbc_pp_err;
        vector<double> runs_pm;
        vector<double> cros_pm;
        vector<double> zdc_pm;
        vector<double> zdc_pm_err;
        vector<double> bbc_pm;
        vector<double> bbc_pm_err;

        TString query =
            "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

        int n_points = t->Draw(query, ppcuts);

        for (std::size_t ip = 0; ip < n_points; ip++) {
          runs_pp.push_back(t->GetV1()[ip]);
          cros_pp.push_back(t->GetV2()[ip]);
          zdc_pp.push_back(t->GetV3()[ip]);
          zdc_pp_err.push_back(t->GetV4()[ip]);
        }

        query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

        int n_points2 = t->Draw(query, ppcuts);

        if (n_points2 != n_points)
          cout << "AAA\n";

        for (std::size_t ip = 0; ip < n_points2; ip++) {
          if (t->GetV1()[ip] != runs_pp[ip])
            cout << "AAA\n";
          if (t->GetV2()[ip] != cros_pp[ip])
            cout << "AAA\n";
          bbc_pp.push_back(t->GetV3()[ip]);
          bbc_pp_err.push_back(t->GetV4()[ip]);
        }

        query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

        n_points = t->Draw(query, pmcuts);

        for (std::size_t ip = 0; ip < n_points; ip++) {
          runs_pm.push_back(t->GetV1()[ip]);
          cros_pm.push_back(t->GetV2()[ip]);
          zdc_pm.push_back(t->GetV3()[ip]);
          zdc_pm_err.push_back(t->GetV4()[ip]);
        }

        query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

        n_points2 = t->Draw(query, pmcuts);

        if (n_points2 != n_points)
          cout << "AAA\n";

        for (std::size_t ip = 0; ip < n_points2; ip++) {
          if (t->GetV1()[ip] != runs_pm[ip])
            cout << "AAA\n";
          if (t->GetV2()[ip] != cros_pm[ip])
            cout << "AAA\n";
          bbc_pm.push_back(t->GetV3()[ip]);
          bbc_pm_err.push_back(t->GetV4()[ip]);
        }

        map<int, int>::iterator rit = m.begin();

        while (rit != m.end()) {
          int currun = rit->first;
          // if (currun > 386946) break;
          TGraphErrors zsame;
          TGraphErrors zoppo;
          TGraphErrors bsame;
          TGraphErrors boppo;

          for (std::size_t ip = 0; ip < runs_pp.size(); ip++) {
            if (runs_pp[ip] == currun) {
              zsame.SetPoint(zsame.GetN(), cros_pp[ip], zdc_pp[ip]);
              zsame.SetPointError(zsame.GetN() - 1, 0, zdc_pp_err[ip]);
              bsame.SetPoint(bsame.GetN(), cros_pp[ip], bbc_pp[ip]);
              bsame.SetPointError(bsame.GetN() - 1, 0, bbc_pp_err[ip]);
            }
          }

          for (std::size_t ip = 0; ip < runs_pm.size(); ip++) {
            if (runs_pm[ip] == currun) {
              zoppo.SetPoint(zoppo.GetN(), cros_pm[ip], zdc_pm[ip]);
              zoppo.SetPointError(zoppo.GetN() - 1, 0, zdc_pm_err[ip]);
              boppo.SetPoint(boppo.GetN(), cros_pm[ip], bbc_pm[ip]);
              boppo.SetPointError(boppo.GetN() - 1, 0, bbc_pm_err[ip]);
            }
          }

          if (zsame.GetN() > 45 && zoppo.GetN() > 45) {

            zsame.Fit("pol0", "q");
            zoppo.Fit("pol0", "q");
            bsame.Fit("pol0", "q");
            boppo.Fit("pol0", "q");

            double zdc_sa = 0.;
            double zdc_sa_err_sq = 0.;
            double zdc_op = 0.;
            double zdc_op_err_sq = 0.;
            double bbc_sa = 0.;
            double bbc_sa_err_sq = 0.;
            double bbc_op = 0.;
            double bbc_op_err_sq = 0.;

            for (std::size_t ip = 0; ip < zsame.GetN(); ip++) {
              zdc_sa += zsame.GetY()[ip];
              zdc_sa_err_sq += pow(zsame.GetEY()[ip], 2);
              bbc_sa += bsame.GetY()[ip];
              bbc_sa_err_sq += pow(bsame.GetEY()[ip], 2);
            }

            for (std::size_t ip = 0; ip < zoppo.GetN(); ip++) {
              zdc_op += zoppo.GetY()[ip];
              zdc_op_err_sq += pow(zoppo.GetEY()[ip], 2);
              bbc_op += boppo.GetY()[ip];
              bbc_op_err_sq += pow(boppo.GetEY()[ip], 2);
            }

            double zdc_sa_err = sqrt(zdc_sa_err_sq);
            double zdc_op_err = sqrt(zdc_op_err_sq);
            double bbc_sa_err = sqrt(bbc_sa_err_sq);
            double bbc_op_err = sqrt(bbc_op_err_sq);

            t->GetEntry(rit->second * 120);
            double pol_blue = t->GetLeaf("gl1p_bpol")->GetValue(0);
            double pol_yell = t->GetLeaf("gl1p_ypol")->GetValue(0);
            double pol_blue_err = t->GetLeaf("gl1p_bpolerr")->GetValue(0);
            double pol_yell_err = t->GetLeaf("gl1p_ypolerr")->GetValue(0);

            double a_ll =
                ((zdc_sa / bbc_sa) - (zdc_op / bbc_op)) /
                (pol_blue * pol_yell * ((zdc_sa / bbc_sa) + (zdc_op / bbc_op)));

            double pol_rel_error_sq = pow(pol_blue_err / pol_blue, 2) +
                                      pow(pol_yell_err / pol_yell, 2);
            double scaler_rel_error_sq =
                pow(zdc_sa_err / zdc_sa, 2) + pow(zdc_op_err / zdc_op, 2) +
                pow(bbc_sa_err / bbc_sa, 2) + pow(bbc_op_err / bbc_op, 2);

            double radicant = pol_rel_error_sq;
            radicant +=
                pow(2 * zdc_sa * zdc_op * bbc_sa * bbc_op, 2) /
                pow(pow(zdc_sa * bbc_op, 2) - pow(zdc_op * bbc_sa, 2), 2) *
                scaler_rel_error_sq;
            double a_ll_err = fabs(a_ll) * sqrt(radicant);

            currresult->SetPoint(currresult->GetN(), currun, a_ll);
            currresult->SetPointError(currresult->GetN() - 1, 0, a_ll_err);
          }
          ++rit;
        }

        currresult->Fit("pol0", "q");
        TF1 *currtf = currresult->GetFunction("pol0");

        std::cout << "currtf->GetParameter(0) " << currtf->GetParameter(0)
                  << std::endl;
        std::cout << "currtf->GetParError(0) " << currtf->GetParError(0)
                  << std::endl;
        // std::cout << "itn " << itn << std::endl;

        result->SetPoint(itn - 1, itn, currtf->GetParameter(0));
        result->SetPointError(itn - 1, 0, currtf->GetParError(0));

        // result->Print();
      }
    }

    return result;
  }

  TGraphErrors *get_cluster_all(const TString yield = "ptbin1arm0", TCut extra_cut = "") {
    TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
                crossing_livetime_cuts + single_arm_qa_cuts +
                north_south_qa_cuts + extra_cut + pileup_cuts + bad_fills +
                fill_cuts + before_spinpats + bad_bunch_shuffling;

    // TCut cuts = final_cuts;

    TCut ppcuts = cuts + "bpat*ypat==1";
    TCut pmcuts = cuts + "bpat*ypat==-1";
    int npp = t->Draw("runnumber", final_cuts + "bpat*ypat==1");
    vector<int> runspp;
    for (std::size_t i = 0; i < npp; i++) {
      runspp.push_back(t->GetV1()[i]);
    }

    int npm = t->Draw("runnumber", final_cuts + "bpat*ypat==-1");
    vector<int> runspm;
    for (std::size_t i = 0; i < npm; i++) {
      runspm.push_back(t->GetV1()[i]);
    }
    // std::cout << "runspp.size() " << runspp.size() << std::endl;
    // std::cout << "runspm.size() " << runspm.size() << std::endl;
    TGraphErrors *result = new TGraphErrors();

    TString bbc_scaler = "fbf_bbcwide * star_clk";
    TString bbc_scaler_err = "fbf_bbcwideerr * star_clk";

    // TString yield = "ptbin";
    // yield += ptbin;
    // // yield += "arm1";
    // yield += "arm0+ptbin";
    // yield += ptbin;
    // yield += "arm1";
    std::cout << "yield " << yield << std::endl;

    vector<double> runs_pp;
    vector<double> cros_pp;
    vector<double> bbc_pp;
    vector<double> bbc_pp_err;
    vector<double> runs_pm;
    vector<double> cros_pm;
    vector<double> bbc_pm;
    vector<double> bbc_pm_err;
    vector<double> yield_pp;
    vector<double> yield_pm;

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    int n_points2 = t->Draw(query, ppcuts);

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      // if (t->GetV1()[ip] != runs_pp[ip]) cout << "AAA\n";
      // if (t->GetV2()[ip] != cros_pp[ip]) cout << "AAA\n";
      runs_pp.push_back(t->GetV1()[ip]);
      cros_pp.push_back(t->GetV2()[ip]);
      bbc_pp.push_back(t->GetV3()[ip]);
      bbc_pp_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    n_points2 = t->Draw(query, pmcuts);

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      // if (t->GetV1()[ip] != runs_pm[ip]) cout << "AAA\n";
      // if (t->GetV2()[ip] != cros_pm[ip]) cout << "AAA\n";
      runs_pm.push_back(t->GetV1()[ip]);
      cros_pm.push_back(t->GetV2()[ip]);
      bbc_pm.push_back(t->GetV3()[ip]);
      bbc_pm_err.push_back(t->GetV4()[ip]);
    }

    query = "runnumber:crossing:" + yield;
    n_points2 = t->Draw(query, ppcuts);

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      yield_pp.push_back(t->GetV3()[ip]);
    }

    n_points2 = t->Draw(query, pmcuts);

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      yield_pm.push_back(t->GetV3()[ip]);
    }

    map<int, int>::iterator rit = m.begin();

    // while (rit != m.end()) {
    while (rit != m.end()) {
    // if (true) {
      int currun = rit->first;
      // if (currun != 398149){
      // 	++rit;
      // 	continue;
      // }

      // if (currun > 386946) break;
      // int currun = 398149;
      int cnpp = 0;
      int cnpm = 0;

      for (std::size_t i = 0; i < npp; i++) {
        if (runspp[i] == currun) {
          cnpp++;
        }
      }

      for (std::size_t i = 0; i < npm; i++) {
        if (runspm[i] == currun) {
          cnpm++;
        }
      }
      // std::cout << "cnpp " << cnpp << std::endl;
      // std::cout << "cnpm " << cnpm << std::endl;

      if (cnpp > 45 && cnpm > 45) {

        double n_pp = 0;
        double n_pm = 0;
	double b_pp = 0;
	double b_pm = 0;
        double R    = 0;

	for (std::size_t ip = 0; ip < runs_pp.size(); ip++) {
            if (runs_pp[ip] == currun // && cros_pp[ip] != 29 && cros_pp[ip] != 30 && cros_pp[ip] != 69 && cros_pp[ip] != 70
		) {
	      b_pp += bbc_pp[ip];
	    }
	}

	for (std::size_t ip = 0; ip < runs_pm.size(); ip++) {
            if (runs_pm[ip] == currun // && cros_pm[ip] != 29 && cros_pm[ip] != 30 && cros_pm[ip] != 69 && cros_pm[ip] != 70
		) {
	      b_pm += bbc_pm[ip];
	    }
	}

	R = b_pp / b_pm;

	for (std::size_t ip = 0; ip < runs_pp.size(); ip++) {
            if (runs_pp[ip] == currun // && cros_pp[ip] != 29 && cros_pp[ip] != 30 && cros_pp[ip] != 69 && cros_pp[ip] != 70
		) {
	      // std::cout << cros_pp[ip] << " yield_pp[ip] " << yield_pp[ip] / bbc_pp[ip] << std::endl;
	      n_pp += yield_pp[ip];
	    }
	}

	// cout << endl;

	for (std::size_t ip = 0; ip < runs_pm.size(); ip++) {
            if (runs_pm[ip] == currun // && cros_pm[ip] != 29 && cros_pm[ip] != 30 && cros_pm[ip] != 69 && cros_pm[ip] != 70
		) {
	      // std::cout << cros_pm[ip] << " yield_pm[ip] " << yield_pm[ip] / bbc_pm[ip] << std::endl;
	      n_pm += yield_pm[ip];
	    }
	}
	// std::cout << "n_pp " << n_pp << std::endl;
	// std::cout << "n_pm " << n_pm << std::endl;

        t->GetEntry(rit->second * 120);
        double pol_blue = t->GetLeaf("gl1p_bpol")->GetValue(0);
        double pol_yell = t->GetLeaf("gl1p_ypol")->GetValue(0);
        double pol_blue_err = t->GetLeaf("gl1p_bpolerr")->GetValue(0);
        double pol_yell_err = t->GetLeaf("gl1p_ypolerr")->GetValue(0);

	// n_pp /= 100;
	// n_pm /= 100;

        double a_ll =
            (n_pp - R * n_pm) / (pol_blue * pol_yell * (n_pp + R * n_pm));
        // double a_ll =
        //     (n_pp - R * n_pm) / ((n_pp + R * n_pm));

        double pol_rel_error_sq =
            pow(pol_blue_err / pol_blue, 2) + pow(pol_yell_err / pol_yell, 2);
        // double scaler_rel_error_sq =
        //     1.0 / n_pp + 1.0 / n_pm + pow((2.616E-5) / R, 2);
        double scaler_rel_error_sq = (1.0 / n_pp) + (1.0 / n_pm);

        // double radicant = pol_rel_error_sq;
	double radicant = 0;
        radicant += pow(2 * R * n_pp * n_pm, 2) /
                    pow(pow(n_pp, 2) - pow(R * n_pm, 2), 2) *
                    scaler_rel_error_sq;
        double a_ll_err = fabs(a_ll) * sqrt(radicant);
	double asym_error = (1.0/(pol_blue * pol_yell)) * 2 * R * n_pp * n_pm / (pow((n_pp + R * n_pm),2));
	asym_error = asym_error*sqrt((1.0 / n_pp) + (1.0 / n_pm));
	// double asym_error = 2 * R * n_pp * n_pm / (pow((n_pp + R * n_pm),2));
	// asym_error = asym_error*sqrt((1.0 / n_pp) + (1.0 / n_pm));

        if (false) {
	  std::cout << "n_pp " << n_pp << std::endl;
	  std::cout << "n_pm " << n_pm << std::endl;
	  std::cout << "R " << R << std::endl;
	  // std::cout << "pol_blue " << pol_blue << std::endl;
	  // std::cout << "pol_yell " << pol_yell << std::endl;
	  std::cout << "scaler_rel_error_sq " << scaler_rel_error_sq << std::endl;
	  std::cout << "sqrt(radicant) " << sqrt(radicant) << std::endl;
	  std::cout << "a_ll " << a_ll << std::endl;
	  std::cout << "a_ll_err " << a_ll_err << std::endl;
	}

        if (a_ll > -1 && a_ll < 1 && a_ll_err > 0 && a_ll_err < 1) {
          result->SetPoint(result->GetN(), currun, a_ll);
          // result->SetPointError(result->GetN() - 1, 0, a_ll_err);
          result->SetPointError(result->GetN() - 1, 0, asym_error);
        }
      }
      ++rit;
    }

    return result;
  }

  TGraphErrors *get_cluster_al(const TString yield = "ptbin1arm0",
                               const bool is_blue, TCut extra_cut = "") {
    TCut cuts = final_cuts + extra_cut;

    TString pat = (is_blue ? "bpat" : "ypat");

    int npp = t->Draw("runnumber", final_cuts + "bpat*ypat==1");
    vector<int> runspp;
    for (std::size_t i = 0; i < npp; i++) {
      runspp.push_back(t->GetV1()[i]);
    }

    int npm = t->Draw("runnumber", final_cuts + "bpat*ypat==-1");
    vector<int> runspm;

    for (std::size_t i = 0; i < npm; i++) {
      runspm.push_back(t->GetV1()[i]);
    }

    TCut pcuts = cuts && TString(pat + "==1");
    TCut mcuts = cuts && TString(pat + "==-1");

    TGraphErrors *result = new TGraphErrors();

    vector<double> runs_p;
    vector<double> cros_p;
    vector<double> yld_p;
    vector<double> bbc_p;
    vector<double> runs_m;
    vector<double> cros_m;
    vector<double> yld_m;
    vector<double> bbc_m;

    TString query = "runnumber:crossing:" + yield;

    int n_points = t->Draw(query, pcuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      runs_p.push_back(t->GetV1()[ip]);
      cros_p.push_back(t->GetV2()[ip]);
      yld_p.push_back(t->GetV3()[ip]);
    }

    TString bbc_scaler = "fbf_bbcwide * star_clk";

    query = "runnumber:crossing:" + bbc_scaler;

    int n_points2 = t->Draw(query, pcuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      if (t->GetV1()[ip] != runs_p[ip])
        cout << "AAA\n";
      if (t->GetV2()[ip] != cros_p[ip])
        cout << "AAA\n";
      bbc_p.push_back(t->GetV3()[ip]);
    }

    query = "runnumber:crossing:" + yield;

    n_points = t->Draw(query, mcuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      runs_m.push_back(t->GetV1()[ip]);
      cros_m.push_back(t->GetV2()[ip]);
      yld_m.push_back(t->GetV3()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler;

    n_points2 = t->Draw(query, mcuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      if (t->GetV1()[ip] != runs_m[ip])
        cout << "AAA\n";
      if (t->GetV2()[ip] != cros_m[ip])
        cout << "AAA\n";
      bbc_m.push_back(t->GetV3()[ip]);
    }

    // std::size_t index_pp = 0;
    // std::size_t index_m = 0;

    map<int, int>::iterator rit = m.begin();

    while (rit != m.end()) {

      int currun = rit->first;
      // if (currun != 398149){
      // 	++rit;
      // 	continue;
      // }

      int cnpp = 0;
      int cnpm = 0;

      for (std::size_t i = 0; i < npp; i++) {
        if (runspp[i] == currun) {
          cnpp++;
        }
      }

      for (std::size_t i = 0; i < npm; i++) {
        if (runspm[i] == currun) {
          cnpm++;
        }
      }

      if (cnpp > 45 && cnpm > 45) {

        double n_p = 0;
        double n_m = 0;
        double b_p = 0;
        double b_m = 0;

        for (std::size_t ip = 0; ip < runs_p.size(); ip++) {
          if (runs_p[ip] == currun) {
            n_p += yld_p[ip];
            b_p += bbc_p[ip];
          }
        }

        for (std::size_t ip = 0; ip < runs_m.size(); ip++) {
          if (runs_m[ip] == currun) {
            n_m += yld_m[ip];
            b_m += bbc_m[ip];
          }
        }
	
	if (true){
	  TString beam = is_blue? "gl1p_bpol" : "gl1p_ypol";
	
	  double r = b_p / b_m;

	  t->GetEntry(rit->second * 120);
	  double pol = t->GetLeaf(beam)->GetValue(0);

	  // std::cout << "n_p " << n_p << std::endl;
	  // std::cout << "n_m " << n_m << std::endl;
	  // std::cout << "b_p " << b_p << std::endl;
	  // std::cout << "b_m " << b_m << std::endl;
	  // std::cout << "r " << r << std::endl;
	  // std::cout << "pol " << pol << std::endl;

	  double a_l = (n_p - r * n_m) / (pol * (n_p + r * n_m));

	  double yield_error_sq = (1. / n_p) + (1. / n_m);

	  double radicant = 0;
	  radicant += pow(2 * r * n_p * n_m, 2) /
	    pow(pow(n_p, 2) - pow(r * n_m, 2), 2) * yield_error_sq;
	  double a_l_err = fabs(a_l) * sqrt(radicant);

	  if (a_l > -1 && a_l < 1 && a_l_err > 0 && a_l_err < 1) {
	    // std::cout << "a_l " << a_l << std::endl;
	    // std::cout << "a_l_err " << a_l_err << std::endl;
	    result->SetPoint(result->GetN(), currun, a_l);
	    result->SetPointError(result->GetN() - 1, 0, a_l_err);
	  }
	}
      }
      ++rit;
    }
    return result;
  }

  TGraphErrors *
  shuffle_run(const int runnumber, const int mode = 0,
	      const double err_factor = 1.) {
    TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
      crossing_livetime_cuts + single_arm_qa_cuts +
      north_south_qa_cuts + pileup_cuts + bad_fills + fill_cuts +
      before_spinpats;

    gROOT->ProcessLine(".L "
                       "/phenix/spin2/pmontu/offline/analysis/pmontu/"
                       "relative_luminosity/macros/shuffle.C+");

    TString runcut = "runnumber==";
    runcut += runnumber;

    cuts += runcut;

    // TH1D *result = new TH1D("result", "bunch shuffling", 101, -6., 6.);

    TString gname = "gbsh";
    gname += runnumber;
    gname += "_";
    gname += mode == 0 ? 0 : mode;

    TString gtitle = "Bunch Shuffling ";
    gtitle += runnumber;
    gtitle += mode == 0 ? " - Random Helicities" : " Restricted Shuffle, ";
    gtitle += mode == 0 ? "" : mode;
    gtitle += mode == 0 ? "" : " group(s)";

    TGraphErrors *result = new TGraphErrors();
    result->SetName(gname);
    result->SetTitle(gtitle);
    result->SetMarkerStyle(21);
    result->SetMarkerSize(0.5);
    result->SetMarkerColor(kBlack);
    result->SetLineColor(kBlack);

    result->GetXaxis()->SetTitle("Iteration");
    result->GetYaxis()->SetTitle("A_LL");

    TString bbc_scaler = "fbf_bbcwide * star_clk";
    TString zdc_scaler = "fbf_zdcwide * star_clk";

    TString bbc_scaler_err = "fbf_bbcwideerr * star_clk";
    TString zdc_scaler_err = "fbf_zdcwideerr * star_clk * ";

    zdc_scaler_err += err_factor;

    // std::cout << "bbc_scaler_err " << bbc_scaler_err << std::endl;

    vector<double> zdc;
    vector<double> zdc_err;
    vector<double> bbc;
    vector<double> bbc_err;
    vector<double> helicity;
    vector<double> crossing;

    TString query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

    int n_points = t->Draw(query, cuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      zdc.push_back(t->GetV3()[ip]);
      zdc_err.push_back(t->GetV4()[ip]);
      crossing.push_back(t->GetV2()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    int n_points2 = t->Draw(query, cuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      bbc.push_back(t->GetV3()[ip]);
      bbc_err.push_back(t->GetV4()[ip]);
    }

    n_points2 = t->Draw("runnumber:crossing:bpat*ypat", cuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      helicity.push_back(t->GetV3()[ip]);
    }

    vector<double> zdc_pp;
    vector<double> zdc_pp_err;
    vector<double> bbc_pp;
    vector<double> bbc_pp_err;
    vector<double> zdc_pm;
    vector<double> zdc_pm_err;
    vector<double> bbc_pm;
    vector<double> bbc_pm_err;

    shuffle sh = shuffle(zdc, bbc, zdc_err, bbc_err, helicity, crossing, mode);

    for (std::size_t i = 0; i < 10000; i++) {

      if (i % 1000 == 0)
        cout << i / 100 << "%" << endl;

      sh.rnd_shuffle(zdc_pp, bbc_pp, zdc_pp_err, bbc_pp_err, zdc_pm, bbc_pm,
                     zdc_pm_err, bbc_pm_err);

      // std::cout << "zdc.size() " << zdc.size() << std::endl;
      // std::cout << "zdc_pp.size() " << zdc_pp.size() << std::endl;
      // std::cout << "zdc_pm.size() " << zdc_pm.size() << std::endl;

      TGraphErrors zsame;
      TGraphErrors zoppo;
      TGraphErrors bsame;
      TGraphErrors boppo;

      for (std::size_t ip = 0; ip < zdc_pp.size(); ip++) {
        zsame.SetPoint(zsame.GetN(), zsame.GetN(), zdc_pp[ip]);
        zsame.SetPointError(zsame.GetN() - 1, 0, zdc_pp_err[ip]);
        bsame.SetPoint(bsame.GetN(), bsame.GetN(), bbc_pp[ip]);
        bsame.SetPointError(bsame.GetN() - 1, 0, bbc_pp_err[ip]);
      }

      for (std::size_t ip = 0; ip < zdc_pm.size(); ip++) {
        zoppo.SetPoint(zoppo.GetN(), zoppo.GetN(), zdc_pm[ip]);
        zoppo.SetPointError(zoppo.GetN() - 1, 0, zdc_pm_err[ip]);
        boppo.SetPoint(boppo.GetN(), boppo.GetN(), bbc_pm[ip]);
        boppo.SetPointError(boppo.GetN() - 1, 0, bbc_pm_err[ip]);
      }

      if (true) {
        // zsame.Fit("pol0", "q");
        // zoppo.Fit("pol0", "q");
        // bsame.Fit("pol0", "q");
        // boppo.Fit("pol0", "q");

        double zdc_sa = 0.;
        double zdc_sa_err_sq = 0.;
        double zdc_op = 0.;
        double zdc_op_err_sq = 0.;
        double bbc_sa = 0.;
        double bbc_sa_err_sq = 0.;
        double bbc_op = 0.;
        double bbc_op_err_sq = 0.;

        for (std::size_t ip = 0; ip < zsame.GetN(); ip++) {
          zdc_sa += zsame.GetY()[ip];
          zdc_sa_err_sq += pow(zsame.GetEY()[ip], 2);
          bbc_sa += bsame.GetY()[ip];
          bbc_sa_err_sq += pow(bsame.GetEY()[ip], 2);
        }

        for (std::size_t ip = 0; ip < zoppo.GetN(); ip++) {
          zdc_op += zoppo.GetY()[ip];
          zdc_op_err_sq += pow(zoppo.GetEY()[ip], 2);
          bbc_op += boppo.GetY()[ip];
          bbc_op_err_sq += pow(boppo.GetEY()[ip], 2);
        }

        double zdc_sa_err = sqrt(zdc_sa_err_sq);
        double zdc_op_err = sqrt(zdc_op_err_sq);
        double bbc_sa_err = sqrt(bbc_sa_err_sq);
        double bbc_op_err = sqrt(bbc_op_err_sq);

        t->GetEntry(get_run_index(runnumber) * 120);
        double pol_blue = t->GetLeaf("gl1p_bpol")->GetValue(0);
        double pol_yell = t->GetLeaf("gl1p_ypol")->GetValue(0);
        double pol_blue_err = t->GetLeaf("gl1p_bpolerr")->GetValue(0);
        double pol_yell_err = t->GetLeaf("gl1p_ypolerr")->GetValue(0);

        double a_ll =
	  ((zdc_sa / bbc_sa) - (zdc_op / bbc_op)) /
	  (pol_blue * pol_yell * ((zdc_sa / bbc_sa) + (zdc_op / bbc_op)));

        double pol_rel_error_sq =
	  pow(pol_blue_err / pol_blue, 2) + pow(pol_yell_err / pol_yell, 2);
        double scaler_rel_error_sq =
	  pow(zdc_sa_err / zdc_sa, 2) + pow(zdc_op_err / zdc_op, 2) +
	  pow(bbc_sa_err / bbc_sa, 2) + pow(bbc_op_err / bbc_op, 2);

        double radicant = pol_rel_error_sq;
        radicant += pow(2 * zdc_sa * zdc_op * bbc_sa * bbc_op, 2) /
	  pow(pow(zdc_sa * bbc_op, 2) - pow(zdc_op * bbc_sa, 2), 2) *
	  scaler_rel_error_sq;
        double a_ll_err = fabs(a_ll) * sqrt(radicant);

        if (false) {
          std::cout << "zdc_sa " << zdc_sa << std::endl;
          std::cout << "zdc_sa_err " << zdc_sa_err << std::endl;
          std::cout << "zdc_op " << zdc_op << std::endl;
          std::cout << "zdc_op_err " << zdc_op_err << std::endl;
          std::cout << "bbc_sa " << bbc_sa << std::endl;
          std::cout << "bbc_sa_err " << bbc_sa_err << std::endl;
          std::cout << "bbc_op " << bbc_op << std::endl;
          std::cout << "bbc_op_err " << bbc_op_err << std::endl;
          std::cout << "pol_blue " << pol_blue << std::endl;
          std::cout << "pol_yell " << pol_yell << std::endl;
          std::cout << "pol_blue_err " << pol_blue_err << std::endl;
          std::cout << "pol_yell_err " << pol_yell_err << std::endl;
          std::cout << "a_ll " << a_ll << std::endl;
          std::cout << "pol_rel_error_sq " << pol_rel_error_sq << std::endl;
          std::cout << "scaler_rel_error_sq " << scaler_rel_error_sq
                    << std::endl;
          std::cout << "radicant " << radicant << std::endl;
          std::cout << "a_ll_err " << a_ll_err << std::endl;
        }
        if (true) {
          // std::cout << "(a_ll / a_ll_err) " << (a_ll / a_ll_err) <<
          // std::endl;

          result->SetPoint(result->GetN(), result->GetN(), a_ll);
          result->SetPointError(result->GetN() - 1, 0, a_ll_err);
        }
        zdc_pp.clear();
        bbc_pp.clear();
        zdc_pp_err.clear();
        bbc_pp_err.clear();
        zdc_pm.clear();
        bbc_pm.clear();
        zdc_pm_err.clear();
        bbc_pm_err.clear();
      }
    }
    return result;
  }

  TGraphErrors *al_shuffle_run(const int runnumber, const int mode = 0,
                               const bool is_blue = true,
                               const double err_factor = 1.) {

    TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
                crossing_livetime_cuts + single_arm_qa_cuts +
                north_south_qa_cuts + pileup_cuts + bad_fills + fill_cuts +
                before_spinpats;

    gROOT->ProcessLine(".L "
                       "/phenix/spin2/pmontu/offline/analysis/pmontu/"
                       "relative_luminosity/macros/shuffle.C+");

    TString runcut = "runnumber==";
    runcut += runnumber;

    cuts += runcut;

    // TH1D *result = new TH1D("result", "bunch shuffling", 101, -6., 6.);

    TString gname = "g_";
    gname += is_blue ? "b" : "y";
    gname += "_al_bsh";
    gname += runnumber;
    gname += "_";
    gname += mode == 0 ? 0 : mode;

    TString gtitle = is_blue ? "Blue" : "Yellow";
    gtitle += " A_L Bunch Shuffling ";
    gtitle += runnumber;
    gtitle += mode == 0 ? " - Random Helicities" : " Restricted Shuffle, ";
    gtitle += mode == 0 ? "" : mode;
    gtitle += mode == 0 ? "" : " group(s)";

    TGraphErrors *result = new TGraphErrors();
    result->SetName(gname);
    result->SetTitle(gtitle);
    result->SetMarkerStyle(21);
    result->SetMarkerSize(0.5);
    result->SetMarkerColor(kBlack);
    result->SetLineColor(kBlack);

    result->GetXaxis()->SetTitle("Iteration");
    result->GetYaxis()->SetTitle("A_L");

    TString bbc_scaler = "fbf_bbcwide * star_clk";
    TString zdc_scaler = "fbf_zdcwide * star_clk";

    TString bbc_scaler_err = "fbf_bbcwideerr * star_clk * ";
    TString zdc_scaler_err = "fbf_zdcwideerr * star_clk * ";

    bbc_scaler_err += err_factor;
    zdc_scaler_err += err_factor;

    // std::cout << "bbc_scaler_err " << bbc_scaler_err << std::endl;

    vector<double> zdc;
    vector<double> zdc_err;
    vector<double> bbc;
    vector<double> bbc_err;
    vector<double> helicity;
    vector<double> crossing;

    TString query = "runnumber:crossing:" + zdc_scaler + ":" + zdc_scaler_err;

    int n_points = t->Draw(query, cuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      zdc.push_back(t->GetV3()[ip]);
      zdc_err.push_back(t->GetV4()[ip]);
      crossing.push_back(t->GetV2()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    int n_points2 = t->Draw(query, cuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      bbc.push_back(t->GetV3()[ip]);
      bbc_err.push_back(t->GetV4()[ip]);
    }

    TString pattern = is_blue ? "bpat" : "ypat";

    n_points2 = t->Draw("runnumber:crossing:" + pattern, cuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      helicity.push_back(t->GetV3()[ip]);
    }

    vector<double> zdc_pp;
    vector<double> zdc_pp_err;
    vector<double> bbc_pp;
    vector<double> bbc_pp_err;
    vector<double> zdc_pm;
    vector<double> zdc_pm_err;
    vector<double> bbc_pm;
    vector<double> bbc_pm_err;

    shuffle sh = shuffle(zdc, bbc, zdc_err, bbc_err, helicity, crossing, mode);

    for (std::size_t i = 0; i < 10000; i++) {

      if (i % 1000 == 0)
        cout << i / 100 << "%" << endl;

      sh.rnd_shuffle(zdc_pp, bbc_pp, zdc_pp_err, bbc_pp_err, zdc_pm, bbc_pm,
                     zdc_pm_err, bbc_pm_err);

      // std::cout << "zdc.size() " << zdc.size() << std::endl;
      // std::cout << "zdc_pp.size() " << zdc_pp.size() << std::endl;
      // std::cout << "zdc_pm.size() " << zdc_pm.size() << std::endl;

      TGraphErrors zsame;
      TGraphErrors zoppo;
      TGraphErrors bsame;
      TGraphErrors boppo;

      for (std::size_t ip = 0; ip < zdc_pp.size(); ip++) {
        zsame.SetPoint(zsame.GetN(), zsame.GetN(), zdc_pp[ip]);
        zsame.SetPointError(zsame.GetN() - 1, 0, zdc_pp_err[ip]);
        bsame.SetPoint(bsame.GetN(), bsame.GetN(), bbc_pp[ip]);
        bsame.SetPointError(bsame.GetN() - 1, 0, bbc_pp_err[ip]);
      }

      for (std::size_t ip = 0; ip < zdc_pm.size(); ip++) {
        zoppo.SetPoint(zoppo.GetN(), zoppo.GetN(), zdc_pm[ip]);
        zoppo.SetPointError(zoppo.GetN() - 1, 0, zdc_pm_err[ip]);
        boppo.SetPoint(boppo.GetN(), boppo.GetN(), bbc_pm[ip]);
        boppo.SetPointError(boppo.GetN() - 1, 0, bbc_pm_err[ip]);
      }

      if (true) {
        // zsame.Fit("pol0", "q");
        // zoppo.Fit("pol0", "q");
        // bsame.Fit("pol0", "q");
        // boppo.Fit("pol0", "q");

        double zdc_sa = 0.;
        double zdc_sa_err_sq = 0.;
        double zdc_op = 0.;
        double zdc_op_err_sq = 0.;
        double bbc_sa = 0.;
        double bbc_sa_err_sq = 0.;
        double bbc_op = 0.;
        double bbc_op_err_sq = 0.;

        for (std::size_t ip = 0; ip < zsame.GetN(); ip++) {
          zdc_sa += zsame.GetY()[ip];
          zdc_sa_err_sq += pow(zsame.GetEY()[ip], 2);
          bbc_sa += bsame.GetY()[ip];
          bbc_sa_err_sq += pow(bsame.GetEY()[ip], 2);
        }

        for (std::size_t ip = 0; ip < zoppo.GetN(); ip++) {
          zdc_op += zoppo.GetY()[ip];
          zdc_op_err_sq += pow(zoppo.GetEY()[ip], 2);
          bbc_op += boppo.GetY()[ip];
          bbc_op_err_sq += pow(boppo.GetEY()[ip], 2);
        }

        double zdc_sa_err = sqrt(zdc_sa_err_sq);
        double zdc_op_err = sqrt(zdc_op_err_sq);
        double bbc_sa_err = sqrt(bbc_sa_err_sq);
        double bbc_op_err = sqrt(bbc_op_err_sq);

        t->GetEntry(get_run_index(runnumber) * 120);
        double pol_blue = t->GetLeaf("gl1p_bpol")->GetValue(0);
        double pol_yell = t->GetLeaf("gl1p_ypol")->GetValue(0);
        double pol_blue_err = t->GetLeaf("gl1p_bpolerr")->GetValue(0);
        double pol_yell_err = t->GetLeaf("gl1p_ypolerr")->GetValue(0);

        double a_ll =
            ((zdc_sa / bbc_sa) - (zdc_op / bbc_op)) /
            (pol_blue * pol_yell * ((zdc_sa / bbc_sa) + (zdc_op / bbc_op)));

        double pol_rel_error_sq =
            pow(pol_blue_err / pol_blue, 2) + pow(pol_yell_err / pol_yell, 2);
        double scaler_rel_error_sq =
            pow(zdc_sa_err / zdc_sa, 2) + pow(zdc_op_err / zdc_op, 2) +
            pow(bbc_sa_err / bbc_sa, 2) + pow(bbc_op_err / bbc_op, 2);

        double radicant = pol_rel_error_sq;
        radicant += pow(2 * zdc_sa * zdc_op * bbc_sa * bbc_op, 2) /
                    pow(pow(zdc_sa * bbc_op, 2) - pow(zdc_op * bbc_sa, 2), 2) *
                    scaler_rel_error_sq;
        double a_ll_err = fabs(a_ll) * sqrt(radicant);

        if (false) {
          std::cout << "zdc_sa " << zdc_sa << std::endl;
          std::cout << "zdc_sa_err " << zdc_sa_err << std::endl;
          std::cout << "zdc_op " << zdc_op << std::endl;
          std::cout << "zdc_op_err " << zdc_op_err << std::endl;
          std::cout << "bbc_sa " << bbc_sa << std::endl;
          std::cout << "bbc_sa_err " << bbc_sa_err << std::endl;
          std::cout << "bbc_op " << bbc_op << std::endl;
          std::cout << "bbc_op_err " << bbc_op_err << std::endl;
          std::cout << "pol_blue " << pol_blue << std::endl;
          std::cout << "pol_yell " << pol_yell << std::endl;
          std::cout << "pol_blue_err " << pol_blue_err << std::endl;
          std::cout << "pol_yell_err " << pol_yell_err << std::endl;
          std::cout << "a_ll " << a_ll << std::endl;
          std::cout << "pol_rel_error_sq " << pol_rel_error_sq << std::endl;
          std::cout << "scaler_rel_error_sq " << scaler_rel_error_sq
                    << std::endl;
          std::cout << "radicant " << radicant << std::endl;
          std::cout << "a_ll_err " << a_ll_err << std::endl;
        }
        if (true) {
          // std::cout << "(a_ll / a_ll_err) " << (a_ll / a_ll_err) <<
          // std::endl;
          // result->Fill(a_ll / a_ll_err);

          result->SetPoint(result->GetN(), result->GetN(), a_ll);
          result->SetPointError(result->GetN() - 1, 0, a_ll_err);
        }
        zdc_pp.clear();
        bbc_pp.clear();
        zdc_pp_err.clear();
        bbc_pp_err.clear();
        zdc_pm.clear();
        bbc_pm.clear();
        zdc_pm_err.clear();
        bbc_pm_err.clear();
      }
    }
    return result;
  }

  // TGraphErrors *get_singles_to_doubles_ratio() {
  //   TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
  //       crossing_livetime_cuts + single_arm_qa_cuts +
  //       north_south_qa_cuts + extra_cut;

  // }

  TGraphErrors *shuffle_cluster_yields(const int runnumber, const int ptbin = 1) {
    
    TCut cuts = final_cuts;

    // gROOT->ProcessLine(".L "
    //                    "/phenix/spin2/pmontu/offline/analysis/pmontu/"
    //                    "relative_luminosity/macros/shuffle.C+");

    TString runcut = "runnumber==";
    runcut += runnumber;

    cuts += runcut;

    // TH1D *result = new TH1D("result", "bunch shuffling", 101, -6., 6.);

    TString gname = "gbsh";
    gname += runnumber;

    TString gtitle = "Bunch Shuffling Run";
    gtitle += runnumber;

    TGraphErrors *result = new TGraphErrors();
    result->SetName(gname);
    result->SetTitle(gtitle);
    result->SetMarkerStyle(21);
    result->SetMarkerSize(0.5);
    result->SetMarkerColor(kBlack);
    result->SetLineColor(kBlack);

    result->GetXaxis()->SetTitle("Iteration");
    result->GetYaxis()->SetTitle("A_LL");

    TString bbc_scaler = "fbf_bbcwide * star_clk";
    TString bbc_scaler_err = "fbf_bbcwideerr * star_clk";
    TString yld_scaler = "ptbin";
    yld_scaler += ptbin;
    yld_scaler += "arm0+ptbin";
    yld_scaler += ptbin;
    yld_scaler += "arm1";

    // std::cout << "bbc_scaler_err " << bbc_scaler_err << std::endl;

    vector<double> yld;
    vector<double> yld_err;
    vector<double> bbc;
    vector<double> bbc_err;
    vector<double> helicity;
    vector<double> crossing;

    TString query = "runnumber:crossing:" + yld_scaler + ":sqrt(" + yld_scaler + ")";
    std::cout << "query " << query << std::endl;

    int n_points = t->Draw(query, cuts);

    for (std::size_t ip = 0; ip < n_points; ip++) {
      yld.push_back(t->GetV3()[ip]);
      yld_err.push_back(t->GetV4()[ip]);
      crossing.push_back(t->GetV2()[ip]);
    }

    query = "runnumber:crossing:" + bbc_scaler + ":" + bbc_scaler_err;

    int n_points2 = t->Draw(query, cuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    for (std::size_t ip = 0; ip < n_points2; ip++) {
      bbc.push_back(t->GetV3()[ip]);
      bbc_err.push_back(t->GetV4()[ip]);
    }

    int fill_n;
    
    n_points2 = t->Draw("runnumber:crossing:bpat*ypat:gl1p_fill", cuts);

    if (n_points2 != n_points)
      cout << "AAA\n";

    fill_n = t->GetV4()[0];
    
    for (std::size_t ip = 0; ip < n_points2; ip++) {
      helicity.push_back(t->GetV3()[ip]);
    }

    vector<double> yld_pp;
    vector<double> yld_pp_err;
    vector<double> bbc_pp;
    vector<double> bbc_pp_err;
    vector<double> yld_pm;
    vector<double> yld_pm_err;
    vector<double> bbc_pm;
    vector<double> bbc_pm_err;

    // Check if random helicity has been generated for this fill
    TString hel_file = "/phenix/spin2/pmontu/offline/analysis/pmontu/"
                       "relative_luminosity/macros/fill_random_helicities/";
    hel_file += fill_n;
    hel_file += ".txt";

    // if (gSystem->AccessPathName(hel_file) == 1) {
    //   ofstream new_hel_file(hel_file);

    //   TRandom3 *r = new TRandom3();
    //   r->SetSeed(0);
    //   for (std::size_t iter = 0; iter < 10000; iter++) {
    // 	for (std::size_t ixing = 0; ixing < 120; ixing++) {
    // 	  new_hel_file << (r->Uniform() > .5 ? -1 : 1) << " ";
    // 	}
    // 	new_hel_file << endl;
    //   }

    //   new_hel_file.close();
    // }

    int rand_helicities[120];

    if (true) {
      TRandom3 *r = new TRandom3();
      r->SetSeed(0);
    }

    ifstream ifile (hel_file);

    for (std::size_t i = 0; i < 10000; i++) {

      if (i % 1000 == 0)
        cout << i / 100 << "%" << endl;

      // for (std::size_t ixing = 0; ixing < 120; ixing++) {
      // 	ifile >> rand_helicities[ixing];
      // 	// std::cout << "rand_helicities[ixing] " << rand_helicities[ixing] << std::endl;
      // }
      for (std::size_t ixing = 0; ixing < 120; ixing++) {
	rand_helicities[ixing] = r->Uniform() > .5 ? -1 : 1;
      }

      for (std::size_t ix = 0; ix < crossing.size(); ix++) {
	if (rand_helicities[int (crossing[ix])] == 1){
	  yld_pp.push_back(yld[ix]);
	  bbc_pp.push_back(bbc[ix]);
	  yld_pp_err.push_back(yld_err[ix]);
	  bbc_pp_err.push_back(bbc_err[ix]);
	} else {
	  yld_pm.push_back(yld[ix]);
	  bbc_pm.push_back(bbc[ix]);
	  yld_pm_err.push_back(yld_err[ix]);
	  bbc_pm_err.push_back(bbc_err[ix]);
	}
      }

      TGraphErrors ysame;
      TGraphErrors yoppo;
      TGraphErrors bsame;
      TGraphErrors boppo;

      for (std::size_t ip = 0; ip < yld_pp.size(); ip++) {
        ysame.SetPoint(ysame.GetN(), ysame.GetN(), yld_pp[ip]);
        ysame.SetPointError(ysame.GetN() - 1, 0, yld_pp_err[ip]);
        bsame.SetPoint(bsame.GetN(), bsame.GetN(), bbc_pp[ip]);
        bsame.SetPointError(bsame.GetN() - 1, 0, bbc_pp_err[ip]);
      }

      for (std::size_t ip = 0; ip < yld_pm.size(); ip++) {
        yoppo.SetPoint(yoppo.GetN(), yoppo.GetN(), yld_pm[ip]);
        yoppo.SetPointError(yoppo.GetN() - 1, 0, yld_pm_err[ip]);
        boppo.SetPoint(boppo.GetN(), boppo.GetN(), bbc_pm[ip]);
        boppo.SetPointError(boppo.GetN() - 1, 0, bbc_pm_err[ip]);
      }

      if (true) {

        double yld_sa = 0.;
        double yld_sa_err_sq = 0.;
        double yld_op = 0.;
        double yld_op_err_sq = 0.;
        double bbc_sa = 0.;
        double bbc_sa_err_sq = 0.;
        double bbc_op = 0.;
        double bbc_op_err_sq = 0.;

        for (std::size_t ip = 0; ip < ysame.GetN(); ip++) {
          yld_sa += ysame.GetY()[ip];
          yld_sa_err_sq += pow(ysame.GetEY()[ip], 2);
          bbc_sa += bsame.GetY()[ip];
          bbc_sa_err_sq += pow(bsame.GetEY()[ip], 2);
        }

        for (std::size_t ip = 0; ip < yoppo.GetN(); ip++) {
          yld_op += yoppo.GetY()[ip];
          yld_op_err_sq += pow(yoppo.GetEY()[ip], 2);
          bbc_op += boppo.GetY()[ip];
          bbc_op_err_sq += pow(boppo.GetEY()[ip], 2);
        }

        double yld_sa_err = sqrt(yld_sa_err_sq);
        double yld_op_err = sqrt(yld_op_err_sq);
        double bbc_sa_err = sqrt(bbc_sa_err_sq);
        double bbc_op_err = sqrt(bbc_op_err_sq);

        t->GetEntry(get_run_index(runnumber) * 120);
        double pol_blue = t->GetLeaf("gl1p_bpol")->GetValue(0);
        double pol_yell = t->GetLeaf("gl1p_ypol")->GetValue(0);
        double pol_blue_err = t->GetLeaf("gl1p_bpolerr")->GetValue(0);
        double pol_yell_err = t->GetLeaf("gl1p_ypolerr")->GetValue(0);

        double a_ll =
	  ((yld_sa / bbc_sa) - (yld_op / bbc_op)) /
	  (pol_blue * pol_yell * ((yld_sa / bbc_sa) + (yld_op / bbc_op)));

        double pol_rel_error_sq =
	  pow(pol_blue_err / pol_blue, 2) + pow(pol_yell_err / pol_yell, 2);
        double scaler_rel_error_sq =
	  pow(yld_sa_err / yld_sa, 2) + pow(yld_op_err / yld_op, 2) +
	  pow(bbc_sa_err / bbc_sa, 2) + pow(bbc_op_err / bbc_op, 2);

        // double radicant = pol_rel_error_sq;
        double radicant = 0;
        radicant += pow(2 * yld_sa * yld_op * bbc_sa * bbc_op, 2) /
	  pow(pow(yld_sa * bbc_op, 2) - pow(yld_op * bbc_sa, 2), 2) *
	  scaler_rel_error_sq;
        double a_ll_err = fabs(a_ll) * sqrt(radicant);

        if (false) {
          std::cout << "yld_sa " << yld_sa << std::endl;
          std::cout << "yld_sa_err " << yld_sa_err << std::endl;
          std::cout << "yld_op " << yld_op << std::endl;
          std::cout << "yld_op_err " << yld_op_err << std::endl;
          std::cout << "bbc_sa " << bbc_sa << std::endl;
          std::cout << "bbc_sa_err " << bbc_sa_err << std::endl;
          std::cout << "bbc_op " << bbc_op << std::endl;
          std::cout << "bbc_op_err " << bbc_op_err << std::endl;
          std::cout << "pol_blue " << pol_blue << std::endl;
          std::cout << "pol_yell " << pol_yell << std::endl;
          std::cout << "pol_blue_err " << pol_blue_err << std::endl;
          std::cout << "pol_yell_err " << pol_yell_err << std::endl;
          std::cout << "a_ll " << a_ll << std::endl;
          std::cout << "pol_rel_error_sq " << pol_rel_error_sq << std::endl;
          std::cout << "scaler_rel_error_sq " << scaler_rel_error_sq
                    << std::endl;
          std::cout << "radicant " << radicant << std::endl;
          std::cout << "a_ll_err " << a_ll_err << std::endl;
        }
        if (true) {
          // std::cout << "(a_ll / a_ll_err) " << (a_ll / a_ll_err) <<
          // std::endl;

          result->SetPoint(result->GetN(), result->GetN(), a_ll);
          result->SetPointError(result->GetN() - 1, 0, a_ll_err);
        }
        yld_pp.clear();
        bbc_pp.clear();
        yld_pp_err.clear();
        bbc_pp_err.clear();
        yld_pm.clear();
        bbc_pm.clear();
        yld_pm_err.clear();
        bbc_pm_err.clear();
      }
    }
    return result;
  }

  int get_run_index(const int run) { return m[run]; }

  // void run_pileup() {
  //   TCut cuts = default_cuts + crossing_qa_cuts + run_livetime_cuts +
  //     crossing_livetime_cuts + single_arm_qa_cuts +
  //     north_south_qa_cuts;
  //   std::cout << "olaKase" << std::endl;
  //   ifstream fr("/phenix/spin2/pmontu/offline/analysis/pmontu/"
  // 		"relative_luminosity/run13_list.txt");

  //   int run, idx == 0;

  //   while (fr >> run && idx < 2) {
  //   std::cout << "run " << run << std::endl;

  //   for (std::size_t ixing = 0; ixing < 2; ixing++) {
  //     TString str = "runnumber == ";
  //     str += run;
  //     str+= " && crossing == ";
  //     str += ixing;
  //     std::cout << (cuts && str).Print() << std::endl;
  //     int n = t->Draw("star_bbcwide", cuts && str, "", 1, idx*120 + ixing);
  //     // std::cout << n << std::endl;
  //     // if (n > 0)
  //     // 	std::cout << str << std::endl;
  //   }

  //   idx++;
  //   }

  // }

private:
  load_map() {
    ifstream fr("/phenix/spin2/pmontu/offline/analysis/pmontu/"
                "relative_luminosity/run13_list.txt");

    int run, index = 0, crossing, entry;
    while (!fr.eof()) {
      fr >> run;
      m[run] = index;
      index++;
    }
    m[398149] = 1086;
  }

  map<int, int> m; // maps from run number to index
};

double true_coincidence_rate(const double south, const double north,
                             const double both) {

  double true_rate =
      log(1.0 - south - north + both) - log(1.0 - south) - log(1.0 - north);

  return true_rate;
}

// // Return true exclusive singles rate from
// // measured inclusive singles rates and true doubles rates
// double true_exclusive_singles_rate(const double singles,
//                                    const double true_coinc_rate) {
//   double true_singles = -log(1.0 - singles) - true_coinc_rate;
//   return true_singles;
// }

// double calc_true_coinc_err(const double clk, const double s, const double n,
// const double both)
// {
//   // cov(N,OR) = <N*OR> - <N><OR> = P_n - P_n*P_or
//   // cov(N,S) = <N*S> - <N><S> = P_both - P_n*P_s
//   Double_t or_cnt = s+n-both;
// //cout << "or_cnt " << or_cnt << endl;
//   Double_t cov_N_OR = n/clk - (n/clk)*(or_cnt/clk);
// //cout << "cov_N_OR " << cov_N_OR << endl;
//   Double_t cov_S_OR = s/clk - (s/clk)*(or_cnt/clk);
// //cout << "cov_S_OR " << cov_S_OR << endl;
//   Double_t cov_N_S = both/clk - (n/clk)*(s/clk);
// //cout << "cov_N_S " << cov_N_S << endl;
//   Double_t sigma_or = (or_cnt/clk)*sqrt(1.0/or_cnt-1.0/clk);    // error is
//   on rate = P_or
// //cout << "sigma_or " << sigma_or << endl;
//   Double_t sigma_n = (n/clk)*sqrt(1.0/n-1.0/clk);    // error is on rate =
//   P_or
// //cout << "sigma_n " << sigma_n << endl;
//   Double_t sigma_s = (s/clk)*sqrt(1.0/s-1.0/clk);    // error is on rate =
//   P_or
// //cout << "sigma_s " << sigma_s << endl;

//   Double_t or_term = -1.0/(1.0-(or_cnt/clk)); // df/dP_or
//   Double_t n_term = 1.0/(1.0-(n/clk));
//   Double_t s_term = 1.0/(1.0-(s/clk));
//   Double_t sigma = or_term*or_term*sigma_or*sigma_or;

//   sigma += n_term*n_term*sigma_n*sigma_n;
//   sigma += s_term*s_term*sigma_s*sigma_s;
//   sigma += 2*or_term*n_term*cov_N_OR*sigma_n*sigma_or;
//   sigma += 2*or_term*s_term*cov_S_OR*sigma_s*sigma_or;
//   sigma += 2*n_term*s_term*cov_N_S*sigma_n*sigma_s;
//   sigma = sqrt(sigma);

//   return sigma;
// }

// double calc_true_exclusive_singles_err(const double clk, const double sgl,
// const double opp, const double both)
// {
//   // cov(N,OR) = <N*OR> - <N><OR> = P_n - P_n*P_or
//   // cov(N,S) = <N*S> - <N><S> = P_both - P_n*P_s
//   double or_cnt = sgl+opp-both;
//   double cov_OPP_OR = (opp/clk - (opp/clk)*(or_cnt/clk));
//   double sigma_or = (or_cnt/clk)*sqrt((1.0/or_cnt-1.0/clk));    // error is
//   on rate = P_or
//   double sigma_opp = (opp/clk)*sqrt((1.0/opp-1.0/clk));

// //  double or_term = -1.0/(1.0-(or_cnt/clk)); // df/dP_or
//   double or_term = 1.0/(1.0-(or_cnt/clk)); // df/dP_or
//   double opp_term = 1.0/(1.0-(opp/clk));

//   double sigma = or_term * or_term * sigma_or * sigma_or;
//   sigma += opp_term * opp_term * sigma_opp * sigma_opp;
//   sigma += 2 * or_term * opp_term * cov_OPP_OR * sigma_or * sigma_opp;
//   sigma = sqrt(sigma);

//   // std::cout << "or_cnt " << or_cnt << std::endl;
//   // std::cout << "cov_OPP_OR " << cov_OPP_OR << std::endl;
//   // std::cout << "sigma_or " << sigma_or << std::endl;
//   // std::cout << "sigma_opp " << sigma_opp << std::endl;
//   // std::cout << "or_term " << or_term << std::endl;
//   // std::cout << "opp_term " << opp_term << std::endl;
//   // std::cout << "sigma " << sigma << std::endl;

//   return sigma;
// }

// very important: these are EXclusive singe-to-double ratios
double FindRate(double r_meas, double kn, double ks) {
  char r_meas_string[200];
  sprintf(r_meas_string, "%f", r_meas);

  TF1 f("Rate Corr", "(-1.0) * [0] + 1.0 - exp(-x * ([1] + 1.0)) - exp(-x*([2] "
                     "+ 1.0))+exp(-x*([1]+[2] + 1.0))",
        0, 1);
  f.SetParameter(0, r_meas);
  f.SetParameter(1, ks);
  f.SetParameter(2, kn);
  ROOT::Math::WrappedTF1 wf1(f);

  ROOT::Math::BrentRootFinder brf;
  brf.SetFunction(wf1, 0, 1);
  brf.Solve();
  if (brf.Status())
    cout << "Error finding root." << endl;
  return brf.Root();
}

double FindRateErr(const double rate, const double nclks) {
  if (nclks <= 0. || rate <= 0.)
    return -999;

  return rate * sqrt((1.0 - rate) / (rate * nclks));
}

double calc_binomial_err(const double rate, const double nclks) {
  if (nclks <= 0. || rate <= 0.)
    return -999;

  return sqrt(rate * (1.0 - rate) / nclks);
}

double calc_coinc_rate_err(double clk, double c, double n, double s) {

  double R_c = c / clk;
  double R_n = n / clk;
  double R_s = s / clk;

  double err_c = calc_binomial_err(R_c, clk);
  double err_n = calc_binomial_err(R_n, clk);
  double err_s = calc_binomial_err(R_s, clk);

  double Partial_c = fabs(1.0 / (1 - R_n - R_s + R_c));
  double Partial_n = fabs(1.0 / (1 - R_n) - Partial_c);
  double Partial_s = fabs(1.0 / (1 - R_s) - Partial_c);

  double var = pow(Partial_c * err_c, 2) + pow(Partial_n * err_n, 2) +
               pow(Partial_s * err_s, 2);

  return sqrt(var);
}

double calc_exclusive_singles_err(const double clk, const double same,
                                  const double oppo, const double coin) {
  double Rsame = same / clk;
  double Roppo = oppo / clk;
  double Rcoin = coin / clk;

  double sameerr = calc_binomial_err(Rsame, clk);
  double oppoerr = calc_binomial_err(Roppo, clk);
  double coinerr = calc_binomial_err(Rcoin, clk);

  double Partial_coin = fabs(1.0 / (1 - Rsame - Roppo + Rcoin));
  double Partial_same = Partial_coin;
  double Partial_oppo = fabs(Partial_coin - 1.0 / (1 - Roppo));

  double var = pow(Partial_coin * coinerr, 2) + pow(Partial_same * sameerr, 2) +
               pow(Partial_oppo * oppoerr, 2);

  return sqrt(var);
}

double calc_k_err(const double clk, const double same, const double oppo,
                  const double coin) {
  // Get rates and errors
  double Rsame = same / clk;
  double Roppo = oppo / clk;
  double Rcoin = coin / clk;

  double sameerr = calc_binomial_err(Rsame, clk);
  double oppoerr = calc_binomial_err(Roppo, clk);
  double coinerr = calc_binomial_err(Rcoin, clk);
  // some terms that show up in error propagation
  double alpha = 1 - Rsame - Roppo + Rcoin;
  double mu = log(1 - Rsame - Roppo + Rcoin) - log(1 - Rsame) - log(1 - Roppo);
  // Partial derivatives
  double Partial_coin = -log(1. - Rsame) / (alpha * pow(mu, 2));
  double Partial_oppo = fabs(((Rsame - Rcoin) * log(1 - Rsame)) /
                             ((1. - Roppo) * alpha * pow(mu, 2)));
  double Partial_same = fabs(1 / (alpha * mu) -
                             (1. / (1. - Rsame) - 1. / alpha) *
                                 (log(1 - Roppo) - log(alpha)) / pow(mu, 2));
  double var = pow(Partial_coin * coinerr, 2) + pow(Partial_same * sameerr, 2) +
               pow(Partial_oppo * oppoerr, 2);
  return sqrt(var);
}

// double FindRateCorr(double r_meas, double kn, double ks) {
//   char r_meas_string[200];
//   sprintf(r_meas_string,"%f",r_meas);

//   // char fstring[200];

//   // sprintf(fstring,"-1*%s + 1 - exp(-x*(%f))
//   -exp(-x*(%f))+exp(-x*(%f+%f-1))",r_meas_string,fks,fkn,fks,fkn);

//   TF1 f("Rate Corr", "(-1.0) * [0] + 1.0 - exp(-x * ([1])) -
//   exp(-x*([2]))+exp(-x*([1]+[2] - 1.0))" , 0, 1);
//   f.SetParameter(0, r_meas);
//   f.SetParameter(1, ks);
//   f.SetParameter(2, kn);
//   ROOT::Math::WrappedTF1 wf1(f);

//   ROOT::Math::BrentRootFinder brf;
//   brf.SetFunction(wf1, 0,1);
//   brf.Solve();
//   if (brf.Status()) cout << "Error finding root." << endl;
//   return brf.Root() / r_meas;
// }

// double FindRateCorrErr(const double nclks, const double nboth)
// {
//   if (nclks == 0. || nboth == 0.)
//     return 0;
//   double rate = nboth / nclks;
//   double err = rate * sqrt((1.0 - rate) / nboth);
//   return err;
// }

// get the rate dependant kn, ks

void xing_by_xing_k(const double coinc_cnt, const double north_cnt,
                    const double south_cnt, const double star_clk,
                    const double coinc_guess, const double guess_err,
                    const string monitor_name, double &true_coinc,
                    double &true_coinc_err, double &true_north,
                    double &true_north_err, double &true_south,
                    double &true_south_err, double &kn, double &kn_err,
                    double &ks, double &ks_err) {

  double c_cnt = coinc_cnt;
  double n_cnt = north_cnt;
  double s_cnt = south_cnt;
  double clock = star_clk;

  bool verb = false;
  // bool verb = true;

  if (verb) {
    std::cout << "c rate " << c_cnt / clock << std::endl;
    std::cout << "n rate " << n_cnt / clock << std::endl;
    std::cout << "s rate " << s_cnt / clock << std::endl;
    std::cout << "clock " << clock << std::endl;
  }

  // first step: correct single arm for pileup
  true_north = -log(1.0 - n_cnt / clock);
  true_south = -log(1.0 - s_cnt / clock);
  true_north_err = calc_binomial_err(n_cnt / clock, clock);
  true_south_err = calc_binomial_err(s_cnt / clock, clock);

  // correct coincidence rate
  true_coinc =
      true_coincidence_rate(s_cnt / clock, n_cnt / clock, c_cnt / clock);
  // true_coinc_err = calc_true_coinc_err(clock, s_cnt, n_cnt, c_cnt);

  true_coinc_err = calc_coinc_rate_err(clock, c_cnt, n_cnt, s_cnt);
  // true_coinc_err = FindRateErr(true_coinc, clock);
  // true_coinc_err = true_coinc * guess_err / coinc_guess;

  // true_coinc = coinc_guess;

  if (verb) {
    std::cout << "true_coinc " << true_coinc << std::endl;
    // std::cout << "corr_coinc_err " << corr_coinc_err << std::endl;
    std::cout << "coinc_guess " << coinc_guess << std::endl;
    std::cout << "guess_err " << guess_err << std::endl;
  }

  // Get the exclusive rates

  double true_excl_north = true_north - true_coinc;
  double true_excl_south = true_south - true_coinc;

  double true_excl_north_err =
      calc_exclusive_singles_err(clock, n_cnt, s_cnt, c_cnt);
  double true_excl_south_err =
      calc_exclusive_singles_err(clock, s_cnt, n_cnt, c_cnt);

  // singles to double ratios:

  kn = true_excl_north / true_coinc;
  ks = true_excl_south / true_coinc;

  kn_err = calc_k_err(clock, n_cnt, s_cnt, c_cnt);
  ks_err = calc_k_err(clock, s_cnt, n_cnt, c_cnt);

  if (verb) {
    std::cout << "true_coinc " << true_coinc << std::endl;
    std::cout << "true_coinc_err " << true_coinc_err << std::endl;
    std::cout << "kn " << kn << std::endl;
    std::cout << "kn_err " << kn_err << std::endl;
    std::cout << "ks " << ks << std::endl;
    std::cout << "ks_err " << ks_err << std::endl;
  }
}

struct grp {
  TString monitor;
  vector<int> run;
  vector<int> crossing;
  vector<double> north;
  vector<double> north_err;
  vector<double> south;
  vector<double> south_err;
  vector<double> uncorr_coinc;
  vector<double> uncorr_coinc_err;
  vector<double> corr_coinc;
  vector<double> corr_coinc_err;
  vector<double> coinc;
  vector<double> coinc_err;
  vector<double> kn;
  vector<double> kn_err;
  vector<double> ks;
  vector<double> ks_err;
  vector<double> rate;
  vector<double> rate_err;
  vector<double> zrate;
  vector<double> zrate_err;
  vector<double> clk;

  void run_iteration(const int itr) {
    vector<double> v(run.size(), 0);
    TGraphErrors *gkn =
        new TGraphErrors(run.size(), &rate[0], &kn[0], &v[0], &kn_err[0]);

    TString name = "g";
    name += monitor;
    name += "knitr";
    name += itr;

    TString title = monitor;
    title += " kn vs Rate Iteration ";
    title += itr;

    gkn->SetName(name);
    gkn->SetTitle(title);
    gkn->SetMarkerStyle(21);
    gkn->SetMarkerSize(0.5);
    gkn->SetMarkerColor(kBlack);
    gkn->SetLineColor(kBlack);

    gkn->GetXaxis()->SetTitle("Rate");
    gkn->GetYaxis()->SetTitle("kn");
    gkn->Fit("pol1", "q");

    TGraphErrors *gks =
        new TGraphErrors(run.size(), &rate[0], &ks[0], &v[0], &ks_err[0]);

    name = "g";
    name += monitor;
    name += "ksitr";
    name += itr;

    title = monitor;
    title += "ks vs Rate Iteration";
    title += itr;

    gks->SetName(name);
    gks->SetTitle(title);
    gks->SetMarkerStyle(21);
    gks->SetMarkerSize(0.5);
    gks->SetMarkerColor(kBlack);
    gks->SetLineColor(kBlack);

    gks->GetXaxis()->SetTitle("Rate");
    gks->GetYaxis()->SetTitle("ks");
    gks->Fit("pol1", "q");

    TF1 *fkn = gkn->GetFunction("pol1");
    TF1 *fks = gks->GetFunction("pol1");
    double currkn = fkn->GetParameter(0);
    double currks = fks->GetParameter(0);

    double currknerr = fkn->GetParError(0);
    double currkserr = fkn->GetParError(0);

    // std::cout << "currkn " << currkn << std::endl;
    // std::cout << "currks " << currks << std::endl;

    for (std::size_t i = 0; i < run.size(); i++) {
      coinc[i] = FindRate(uncorr_coinc[i], currkn, currks);
      kn[i] = (north[i] - coinc[i]) / coinc[i];
      ks[i] = (south[i] - coinc[i]) / coinc[i];
    }
  }

  void save_graphs(const int itr, const TFile *f) {
    f->cd();
    vector<double> v(run.size(), 0);
    TGraphErrors *gkn =
        new TGraphErrors(run.size(), &rate[0], &kn[0], &v[0], &kn_err[0]);

    TString name = "graph";
    name += monitor;
    name += "knitr";
    name += itr;

    TString title = monitor;
    title += " kn vs Rate Iteration ";
    title += itr;

    gkn->SetName(name);
    gkn->SetTitle(title);
    gkn->SetMarkerStyle(21);
    gkn->SetMarkerSize(0.5);
    gkn->SetMarkerColor(kBlack);
    gkn->SetLineColor(kBlack);

    gkn->GetXaxis()->SetTitle("Rate");
    gkn->GetYaxis()->SetTitle("Kn");
    gkn->Fit("pol1", "q");
    gkn->Write();

    TGraphErrors *gks =
        new TGraphErrors(run.size(), &rate[0], &ks[0], &v[0], &ks_err[0]);

    name = "graph";
    name += monitor;
    name += "ksitr";
    name += itr;

    TString title = monitor;
    title += " ks vs Rate Iteration ";
    title += itr;

    gks->SetName(name);
    gks->SetTitle(title);
    gks->SetMarkerStyle(21);
    gks->SetMarkerSize(0.5);
    gks->SetMarkerColor(kBlack);
    gks->SetLineColor(kBlack);

    gks->GetXaxis()->SetTitle("Rate");
    gks->GetYaxis()->SetTitle("Ks");
    gks->Fit("pol1", "q");
    gks->Write();

    // for (std::size_t ip = 0; ip < ndf+2; ip++) {
    //   kn_err[ip] = gkn->GetEY()[ip];
    //   ks_err[ip] = gks->GetEY()[ip];
    // }

    // gks = new TGraphErrors(run.size(), &rate[0], &ks[0], &v[0], &ks_err[0]);
  }
  void get_kn_ks(double &dkn, double &dks, double &dknerr, double &dkserr) {
    vector<double> v(run.size(), 0);
    TGraphErrors *gkn =
        new TGraphErrors(run.size(), &rate[0], &kn[0], &v[0], &kn_err[0]);
    gkn->Fit("pol1");
    TF1 *fn = gkn->GetFunction("pol1");
    int ndf = gkn->GetN() - 2;

    gkn->Fit("pol1");
    fn = gkn->GetFunction("pol1");

    TGraphErrors *gks =
        new TGraphErrors(run.size(), &rate[0], &ks[0], &v[0], &ks_err[0]);
    gks->Fit("pol1");
    TF1 *fs = gks->GetFunction("pol1");
    ndf = gks->GetN() - 2;

    gks->Fit("pol1");
    fs = gks->GetFunction("pol1");

    // helper histograms to estimate RMS

    TH1D *hn = new TH1D(monitor + "hn", monitor + "kN Residuals", 301, -.2, .2);
    TH1D *hs = new TH1D(monitor + "hs", monitor + "kS Residuals", 301, -.2, .2);
    TH1D *whn = new TH1D(monitor + "whn", monitor + "kN Weighted Residuals",
                         301, -.2, .2);
    TH1D *whs = new TH1D(monitor + "whs", monitor + "kS Weighted Residuals",
                         301, -.2, .2);

    for (std::size_t i = 0; i < ndf + 2; i++) {
      double predn = fn->Eval(gkn->GetX()[i]);
      double preds = fs->Eval(gks->GetX()[i]);

      double weightn = pow(gkn->GetEY()[i], -2);
      double weights = pow(gks->GetEY()[i], -2);

      hn->Fill(gkn->GetY()[i] - predn);
      hs->Fill(gks->GetY()[i] - preds);
      whn->Fill(gkn->GetY()[i] - predn, weightn);
      whs->Fill(gks->GetY()[i] - preds, weights);
    }

    dkn = fn->GetParameter(0);
    // dknerr = fn->GetParError(0);
    dknerr = whn->GetRMS();
    dks = fs->GetParameter(0);
    // dkserr = fs->GetParError(0);
    dkserr = whs->GetRMS();

    hn->Write();
    hs->Write();
    whn->Write();
    whs->Write();

    // gks = new TGraphErrors(run.size(), &rate[0], &ks[0], &v[0], &ks_err[0]);
  }
};
