#ifndef BLOCK_HH
#define BLOCK_HH
#include <string>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

// Corsika reader

/*naming conventions of corsika particles:
//
//    1   gamma           24   Omega           64   K* -
//    2   positron        25   anti neutron    65   anti K* 0
//    3   electron        26   anti Lambda     66   electron neutrino
//    4   neutrino        27   anti Sigma -    67   electron anti neutrino
//    5   muon +          28   anti Sigma 0    68   muon neutrino
//    6   muon -          29   anti Sigma +    69   muon anti neutrino
//    7   pion 0          30   anti Xi 0       71   eta-> 2*gam
//    8   pion +          31   anti Xi +       72   eta-> 3*pi0
//    9   pion -          32   anti Omega      73   eta-> pi+ + pi- + pi0
//   10   Kaon 0 long     50   omega           74   eta-> pi+ + pi- + gam
//   11   Kaon            51   rho 0           201   Deuteron
//   12   Kaon -          52   rho +           301   Tritium
//   13   neutron         53   rho -           402   alpha
//   14   proton          54   Delta ++       1206   Carbon
//   15   anti proton     55   Delta +        1407   Nitrogen
//   16   Kaon 0 short    56   Delta 0        1608   Oxygen
//   17   eta (71..74)    57   Delta -        2713   Aluminium
//   18   Lambda          58   anti Delta --  2814   Silicon
//   19   Sigma +         59   anti Delta -   3216   Sulfur
//   20   Sigma 0         60   anti Delta 0   5626   Iron
//   21   Sigma -         61   anti Delta +   9900   Cherenkov photons
//   22   Xi 0            62   K* 0
//   23   Xi -            63   K* +
//  116   D 0            131   tau +           150   anti Xi c -
//  117   D +            132   tau -           151   anti Xi c 0
//  118   anti D -       133   tau neutrino    152   anti Sigma c --
//  119   anti D 0       134   anti tau neutr  153   anti Sigma c -
//  120   D s +          137   Lambda c +      154   anti Sigma c 0
//  121   anti D s -     138   Xi c +          155   anti Xi c prime -
//  122   eta c          139   Xi c 0          156   anti Xi c prime 0
//  123   D*0            140   Sigma c ++      157   anti Omega c 0
//  124   D*+            141   Sigma c +       161   Sigma c * ++
//  125   anti D*-       142   Sigma c 0       162   Sigma c * +
//  126   anti D*0       143   Xi c prime +    163   Sigma c * 0
//  127   D* s +         144   Xi c prime 0    171   anti Sigma c * --
//  128   anti D* s -    145   Omega c 0       172   anti Sigma c * -
//  130   J/psi          149   anti Lambda c-  173   anti Sigma c * 0
----------------------------------------------------------------------- */

/*struct RUNH{
    (1, 'run_header', dtype='S4'),
    (2, 'run_number'),
    (3, 'date'),
    (4, 'version'),
    (5, 'n_observation_levels'),
    (6, 'observation_height', shape=10),
    (16, 'energy_spectrum_slope'),
    (17, 'energy_min'),
    (18, 'energy_max'),
    (19, 'egs4_flag'),
    (20, 'nkg_flag'),
    (21, 'energy_cutoff_hadrons'),
    (22, 'energy_cutoff_muons'),
    (23, 'energy_cutoff_electrons'),
    (24, 'energy_cutoff_photons'),
    (25, 'physical_constants_and_interaction_flags', shape=50),
    (94 + 1, 'cka', shape=40),
    (134 + 1, 'ceta', shape=5),
    (139 + 1, 'cstrba', shape=11),
    (254 + 1, 'aatm', shape=5),
    (259 + 1, 'batm', shape=5),
    (264 + 1, 'catm', shape=5),
    (270, 'nflain'),
    (271, 'nfdif'),
    (272, 'nflpi0_100nflpif'),
    (273, 'nflche_100nfragm'),
}; */

/*struct EVTH{
     (1, 'event_header'),
     (2, 'event_number'),
     (3, 'particle_id'),
     (4, 'total_energy'),
     (5, 'starting_altitude'),
     (6, 'first_target_id'),
     (7, 'first_interaction_height'),
     (8, 'momentum_x'),
     (9, 'momentum_y'),
     (10, 'momentum_minus_z'),
     (11, 'zenith'),
     (12, 'azimuth'),
     (13, 'n_random_sequences'),
     (13 + 1, 'random_seeds'),
     (44, 'run_number'),
     (45, 'date'),
     (46, 'version'),
     (47, 'n_observation_levels'),
     (47 + 1, 'observation_height', shape=10),
     (58, 'energy_spectrum_slope'),
     (59, 'energy_min'),
     (60, 'energy_max'),
     (61, 'energy_cutoff_hadrons'),
     (62, 'energy_cutoff_muons'),
     (63, 'energy_cutoff_electrons'),
     (64, 'energy_cutoff_photons'),
     (65, 'nflain'),
     (66, 'nfdif'),
     (67, 'nflpi0'),
     (68, 'nflpif'),
     (69, 'nflche'),
     (70, 'nfragm'),
     (71, 'earth_magnetic_field_x'),
     (72, 'earth_magnetic_field_z'),
     (73, 'egs4_flag'),
     (74, 'nkg_flag'),
     (75, 'low_energy_hadron_model'),
     (76, 'high_energy_hadron_model'),
     (77, 'cerenkov_flag'),
     (78, 'neutrino_flag'),
     (79, 'curved_flag'),
     (80, 'computer'),
     (81, 'theta_min'),
     (82, 'theta_max'),
     (83, 'phi_min'),
     (84, 'phi_max'),
     (85, 'cherenkov_bunch_size'),
     (86, 'n_cherenkov_detectors_x'),
     (87, 'n_cherenkov_detectors_y'),
     (88, 'cherenkov_detector_grid_spacing_x'),
     (89, 'cherenkov_detector_grid_spacing_y'),
     (90, 'cherenkov_detector_length_x'),
     (91, 'cherenkov_detector_length_y'),
     (92, 'cherenkov_output_flag'),
     (93, 'angle_array_x_magnetic_north'),
     (94, 'additional_muon_information_flag'),
     (95, 'egs4_multpliple_scattering_step_length_factor'),
     (96, 'cherenkov_wavelength_min'),
     (97, 'cherenkov_wavelength_max'),
     (98, 'n_reuse'),
     (98 + 1, 'reuse_x', shape=20),
     (118 + 1, 'reuse_y', shape=20),
     (139, 'sybill_interaction_flag'),
     (140, 'sybill_cross_section_flag'),
     (141, 'qgsjet_interaction_flag'),
     (142, 'qgsjet_cross_section_flag'),
     (143, 'dpmjet_interaction_flag'),
     (144, 'dpmjet_cross_section_flag'),
     (145, 'venus_nexus_epos_cross_section_flag'),
     (146, 'muon_multiple_scattering_flag'),
     (147, 'nkg_radial_distribution_range'),
     (148, 'energy_fraction_if_thinning_level_hadronic'),
     (149, 'energy_fraction_if_thinning_level_em'),
     (150, 'actual_weight_limit_thinning_hadronic'),
     (151, 'actual_weight_limit_thinning_em'),
     (152, 'max_radius_radial_thinning_cutting'),
     (153, 'viewcone_inner_angle'),
     (154, 'viewcone_outer_angle'),
     (155, 'transition_energy_low_high_energy_model'),
};*/

/*struct EVTE{
    EventNumber () const {return fSubBlockData [1];}
    Photons () const {return fSubBlockData [2];}
    Electrons () const {return fSubBlockData [3];}
    Hadrons () const {return fSubBlockData [4];}
    Muons () const {return fSubBlockData [5];}
    Particles () const {return fSubBlockData [6];}
};*/

/*struct RUNE{
    CREAL GetRunNumber () const {return fSubBlockData [1];}
    CREAL GetEventsProcessed () const {return fSubBlockData [2];}
};*/

/*struct DATA{
    int id;
    float px;
    float py;
    float pz;
    float x;
    float y;
    float t;
};*/

class Block{
public:
//    Block(float _radius_cut, float _energy_cut)
//        : radius_cut(_radius_cut), energy_cut(_energy_cut) {}

    static constexpr int block_size = 273;
    float buffer[block_size];
    int total_nevents{};

    void fill_block(std::string header);
    void print_header(std::string path, std::string head);

    // RUN HEADER
    float get_energy_min();
    float get_energy_max();
    float get_nevents();
    float rh_get_element(int index);

    // EVENT HEADER
    float get_primary_id();
    float get_total_energy();
    float get_zenith();
    float get_azimuth();
    float eh_get_element(int index);

    // EVENT END
    float get_photons();
    float get_electrons();
    float get_hadrons();
    float get_muons();
    float get_particles();
    float ee_get_element(int index);

    // RUN END
    float get_run_number();
    float get_n_events();
    float re_get_element(int index);

    // LONG (unused)

    // DATA
      struct Particle {
      Particle() = default;
      Particle(int _id_cor, float _px, float _py, float _pz, float _x, float _y, float _t) { // for raw Corsika file
        //(part.x*part.x+part.y*part.y)<72250000 && part.p>1) Corsika_reader cut
        //energy = std::pow(((_px*_px+_py*_py+_pz*_pz) + gParticleMass[_id_cor/1000]),0.5);
//        if((_x*_x+_y*_y) < radius_cut*radius_cut && (_px*_px+_py*_py+_pz*_pz) + gParticleMass[_id_cor/1000] > energy_cut*energy_cut ){                   // filter cut
//        part.insert(part.end(), {_id_cor/1000, static_cast<int>(std::pow(energy_2*1000, 0.5)),            // transform to geant4 id && to eV
//                                static_cast<int>(_x),static_cast<int>(_y),static_cast<int>(_t)});
             id = _id_cor/1000;
             energy = std::sqrt((_px*_px+_py*_py+_pz*_pz) /*+ gParticleMass[_id_cor/1000]*/)*1000;
             x = _x;
             y = _y;
             t = _t;
//        }
     }
      Particle(float _id_cor, float _en, float _px, float _py, float _pz, float _x, float _y, float _t) { // for compressed file
          id = _id_cor;
          energy = _en;     // MeV
          px = _px;         // MeV
          py = _py;
          pz = _pz;
          x = _x;
          y = _y;
          t = _t;
      }
      int id{};
      int energy{};
      float px{};
      float py{};
      float pz{};
      int x;
      int y;
      int t;
    };
      void print_size(){                                    //NSHOW после фильтров
          std::cout << data.size() << '\n';
      }
      int get_data_size(){
          return data.size();
      }
 
      void print_data(std::string path);
      std::vector<Particle> data;
      std::vector <int> data_block_sizes;
      // other
      float get_mass(int index){
          return gParticleMass[index];
      }
      std::string get_name(int index){
          if (index < 65) {
              return id_s[index];
          }
          else return 0;
      }
      int get_pdg(int index) {
          if (index < 65) {
              return pdg_s[index];
          }
          else return 0;
      }
      void read_file(std::string filename);
      void read_event(std::string filename, int event);
      void read_file_comp(std::string filename);
      void read_event_comp(std::string filename, int event);
      void clear_all();
      int get_nevent_test() {
          return num_of_event;
      }
   /*     void add_params();
        void print_params() {
          for (int i = 0; i < params.size(); ++i) {
              if ((i + 1) % 7 == 0) std::cout << '\n';
              std::cout << params[i] << ' ';
          }
      }*/ 
      void set_core(std::vector<float>& core);
      //void add_params();
      void conv_to_format(std::vector<float> bar, std::string filename, int *Ns, int *N_th);


private:
    static constexpr float radius_cut = 1500.; // cm
    static constexpr float energy_cut = 1.; // MeV
    static constexpr int number_of_det = 16;
    static constexpr int treshold = 6;
    static constexpr int krat = 2;
    static constexpr float mip = 45; // 45 kev
    static constexpr float rp = 1.46;
    float gParticleMass [200] =   //GeV
         {0.e0	,.51099892e-3,.51099892e-3,  0.e0     ,.105658369e0,
         .105658369e0, .1349766e0, .13957018e0,.13957018e0, 0.497648e0 ,//10
         0.493677e0 , 0.493677e0 ,.93956536e0 ,.93827203e0,.93827203e0 ,
         0.497648e0 , 0.54775e0  , 1.115683e0 , 1.18937e0 , 1.192642e0 ,//20
         1.197449e0 , 1.31483e0  , 1.32131e0  , 1.67245e0 ,.93956536e0 ,
         1.115683e0 , 1.18937e0  , 1.192642e0 , 1.197449e0, 1.31483e0  ,//30
         1.32131e0  , 1.67245e0  , 0.e0	  , 0.e0      , 0.e0	   ,
         0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//40
         1.e9	, 580.e0     , 1.e5	  , 0.e0      , 0.e0	   ,
         0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.78259e0  ,//50
         0.7690e0	, 0.7665e0   , 0.7665e0   , 1.2305e0  , 1.2318e0   ,
         1.2331e0	, 1.2344e0   , 1.2309e0   , 1.2323e0  , 1.2336e0   ,//60
         1.2349e0	, 0.89610e0  , 0.89166e0  , 0.89166e0 , 0.89610e0  ,
         0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//70
         0.54775e0  , 0.54775e0  , 0.54775e0  , 0.54775e0 ,.105658369e0,
         .105658369e0 , 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//80
         0.e0	, 0.e0       , 0.e0	  , 0.e0      ,.105658369e0,
         .105658369e0, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//90
         0.e0	, 0.e0       , 0.e0	  , 0.e0      ,.105658369e0,
         .105658369e0, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//100
         0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,
         0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,//110
         0.e0	, 0.e0       , 0.e0	  , 0.e0      , 0.e0	   ,
         1.8645e0	, 1.8697e0   , 1.8697e0   , 1.8645e0   , 1.9682e0   ,//120
         1.9682e0	, 2.9804e0   , 2.0067e0   , 2.0100e0   , 2.0100e0   ,
         2.0067e0	, 2.1121e0   , 2.1121e0   , 0.0e0      , 3.096916e0 ,//130
         1.77699e0  , 1.77699e0  , 0.e0	  , 0.e0       , 0.e0	    ,
         0.e0	, 2.28646e0  , 2.4679e0   , 2.4710e0   , 2.45402e0  ,//140
         2.4529e0	, 2.45376e0  , 2.5757e0   , 2.5780e0   , 2.6975e0   ,
         0.e0	, 0.e0       , 0.e0	  , 2.28646e0  , 2.4679e0   ,//150
         2.4710e0	, 2.45402e0  , 2.4529e0   , 2.45376e0  , 2.5757e0   ,
         2.5780e0	, 2.6975e0   , 0.e0	  , 0.e0       , 0.e0	    ,//160
         2.5184e0	, 2.5175e0   , 2.5180e0   , 0.e0       , 0.e0	    ,
         0.e0	, 0.e0       , 0.e0	  , 0.e0       , 0.e0	    ,//170
         2.5184e0	, 2.5175e0   , 2.5180e0   , 0.e0       , 0.e0	    ,
         5.27958e0  , 5.27925e0  , 5.27925e0  , 5.27958e0  , 5.36677e0  ,//180
         5.36677e0  , 6.277e0    , 6.277e0    , 5.6194e0   , 5.8155e0   ,
         5.8113e0	, 5.788e0    , 5.7911e0   , 6.071e0    , 5.6194e0   ,//190
         5.8155e0	, 5.8113e0   , 5.788e0    , 5.7911e0   , 6.071e0    ,
         0.e0	, 0.e0       , 0.e0	  , 0.e0       , 0.e0	    };
    std::string id_s[202] =
        {"", "gamma", "e+", "e-", "", "mu+", "mu-", "pi0", "pi+", "pi-", // 10
         "K0_L", "K+", "K-", "n", "p", "p_bar", "K0_S", "eta", "Lambda", "Sigma+", //20
         "Sigma 0", "Sigma-", "Xi0", "Xi-", "Omega-", "n_bar", "Lambda_bar", "Sigma-_bar", "Sigma0_bar", "Sigma+_bar", "Xi0_bar", //30
         "Xi+_bar","Omega+_bar", "", "", "", "", "", "", "", "", //40
         "", "", "", "", "", "", "","","", "omega",//50
         "rho0", "rho+","rho-", "Delta++", "Delta+", "Delta0", "Delta-", "Delta--_bar", "Delta-_bar", "Delta0_bar", // 60.
         "Delta+_bar","K*0","K*+", "K*-", "K*0_bar","nu_e",
         "nu_e_bar", "nu_mu", "nu_mu_bar", "", // 70
         "eta->2gamma","eta->3pi0", "eta->pi+pi-pi0", "eta->pi+pi-gamma", "mu+ prod. info", "mu- prod. info", 
         "","", "","", //80
         "", "","","", "dky mu+ start info", "dky mu- start info", "", "","","", // 90.
         "","","", "","dky mu+ end info","dky mu- end info", "","","", "", //100
         "","", "","","", "","","", "","",//110
        "", "", "","", "","D0","D+", "D-_bar", "D0_bar", "D_a+", // 120.

         "D_a-_bar","eta_c", "D*0", "D*+","D*-_bar", " anti D * 0", " D * s +"," anti D * s -"," ",
         " J/psi"," tau +"," tau -", " tau neutrino"," anti tau neutrino "," ", " "," Lambda c +"," Xi c + ",
         " Xi c 0 "," Sigma c ++"," Sigma c +", " Sigma c 0"," Xi c prime +"," Xi c prime 0", " Omega c 0"," "," ",
         " "," anti Lambda c -"," anti Xi c - ", // 121 to 150.
         " anti Xi c 0 "," anti Sigma c --"," anti Sigma c - ", " anti Sigma c 0 "," anti Xi c prime - "," anti Xi c prime 0 ",
         " anti Omega c 0 "," "," ", " "," Sigma c * ++"," Sigma c * + ", " Sigma c * 0 "," "," ", " "," "," ", " "," "," anti Sigma c * -- ",
         " anti Sigma c * -"," anti Sigma c * 0"," ", " "," "," ", " "," "," ", // 151 to 180.
         " "," "," ", " "," "," ", " "," "," ", " "," "," ", " "," "," ", " "," "," ", " Cherenkov photon"," "," "};
    int pdg_s[201] = {0, 22, -11, 11, 0, -13,
    13, 111, 211, -211, 130,              // 10
    321, -321, 2112, 2212, -2212,
    310, 221, 3122, 3222, 3212,        // 20
    3112, 3322, 3312, 3334, -2112,
    -3122, -3112, -3212, -3222, -3322, // 30
    -3312, -3334, 0, 0, 0,
    0, 0, 0, 0, 0,                     // 40
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 223,                     // 50
    113, 213, -213, 2224, 2214,
    2114, 1114, -2224, -2214, -2114,     // 60
    -2214, 313, 323, -323, -313,
    12, -12, 14, -14, 0,                 // 70
    0, 0, 0, 0, -13,
    13, 0, 0, 0, 0,                        // 80
    0, 0, 0, 0, -13,
    13, 0, 0, 0, 0,                        // 90
    0, 0, 0, 0, -13,
    13, 0, 0, 0, 0,                        // 100
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,                        // 110
    0, 0, 0, 0, 0,
    421, 411, -411, -421, 431,           // 120
    -431, 441, 423, 413, -413,
    -423, 433, -433, 0, 443,              // 130
    -15, 15, 16, -16, 0,
    0, 4122, 4232, 4132, 4222,           // 140
    4212, 4112, 4322, 4312, 4332,
    0,0,0, -4122, -4323,                 // 150
    -4132, -4222, -4212, -4112, -4322,
    -4312, -4332, 0, 0, 0,                // 160
    4224, 4214, 4114, 0, 0,
    0, 0, 0, 0, 0,                        // 170
    -4224, -4214, -4114, 0, 0,
    0, 0, 0, 0, 0,                        // 180
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,                        // 190
    0, 0, 0, 0, 0,
    0, 0, 0, 0, 0                        // 200
    };

    int num_of_event{};

    std::vector<float> run_header;
    std::vector<float> event_header;
    std::vector<float> event_end;
    std::vector<float> run_end;
    std::vector<float> data_block;
    std::vector<double> params;
};

#endif // BLOCK_HH
