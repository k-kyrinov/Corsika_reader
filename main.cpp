//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <sstream>
//#include <string>
//#include <cstring>
#include "block.h"

//struct Particle {
//Particle() = default;
//Particle(int _id_cor, float _px, float _py, float _pz, float _x, float _y, float _t) {
//  //(part.x*part.x+part.y*part.y)<72250000 && part.p>1) Corsika_reader cut
//  //energy = std::pow(((_px*_px+_py*_py+_pz*_pz) + gParticleMass[_id_cor/1000]),0.5);
////        if((_x*_x+_y*_y) < radius_cut*radius_cut && (_px*_px+_py*_py+_pz*_pz) + gParticleMass[_id_cor/1000] > energy_cut*energy_cut ){                   // filter cut
////        part.insert(part.end(), {_id_cor/1000, static_cast<int>(std::pow(energy_2*1000, 0.5)),            // transform to geant4 id && to eV
////                                static_cast<int>(_x),static_cast<int>(_y),static_cast<int>(_t)});
//       id = _id_cor/1000;
//       energy = std::sqrt((_px*_px+_py*_py+_pz*_pz) /*+ gParticleMass[static_cast<int>(_id_cor/1000)]*/)*1000;
//       x = _x;
//       y = _y;
//       t = _t;
////        }
//}
//int id{};
//int energy{}; // MeV
//int x;        // cm
//int y;
//int t;        // ns
//};


int main(/*int argc, char *argv[]*/)
{
    Block bl;
/*    std::ifstream file;  //("F://1.New//Corsika//p500TeVtest_1");
    std::ofstream output ("F://1.New//Corsika//test_out.txt");
    float sizeBuff;
    float buffer[273];
    file.open("F://1.New//Corsika//p500TeVtest_1", std::fstream::binary | std::fstream::out); //p500TeVtest_1
    int num_of_event = 0;
    int event = 1;
    std::vector<float> data_info;
    std::vector<Particle> data;
    int block_size = 273;
    float radius_cut = 1500.; // cm       filter cut
    float energy_cut = 1.; // MeV
    int j = 0;
    if(file){
        file.seekg(0);
        file.read((char *)&sizeBuff, 4);
        while(!file.eof() ){
        file.read((char *)&buffer, sizeof (buffer));
        for(int i = 0; i < 273; ++i) {
            if(i==0){
                std::string oString;
                oString = reinterpret_cast<char*>(&buffer[0]);
                if(oString.compare("EVTH")==0){
                    data_info.insert(data_info.end(), {buffer[2], buffer[3], buffer[10], buffer[11]});
                    num_of_event+=1;
                }
                if(event==num_of_event) {
                    if(oString.compare("RUNH")!=0 && oString.compare("EVTH")!=0 && oString.compare("RUNE")!=0 ){
                    for(int i = 0; i <= block_size-7; i+=7){
                        if(buffer[i]!=0 && (buffer[i+4]*buffer[i+4]+buffer[i+5]*buffer[i+5]) < radius_cut*radius_cut &&
                        std::sqrt((buffer[i+1]*buffer[i+1]+buffer[i+2]*buffer[i+2]+buffer[i+3]*buffer[i+3]))*1000. > energy_cut){
                        Particle p(buffer[i], buffer[i+1], buffer[i+2], buffer[i+3], buffer[i+4], buffer[i+5], buffer[i+6]);
                        data.push_back(p);
                        }
                    }
                }
             }
                if(num_of_event>event || oString.compare("EVTE")==0 ) break;
            }
        }

        j+=1;
        if(j%21==0){
            for(int i = 0; i < 2; ++i){
            file.read((char *)&sizeBuff, 4);
            }
        }
     }

        std::ofstream output("F://1.New//Corsika//output_data.txt");
        output << data_info[0] << ' ' << data_info[1] << ' ' << data_info[2] << ' ' << data_info[3] << '\n';
        for(int i = 0; i < data.size(); ++i){
            output << data[i].id << ' ' << data[i].energy << ' '<< data[i].x << ' '<< data[i].y << ' '<< data[i].t << '\n';
        }

    }

    std::cout << data.size() << ' ' << data_info.size() <<'\n'; */


//    bl.read_file("F://1.New//Corsika//p500TeVtest_1");
//    //bl.read_event("F://1.New//Corsika//p500TeVtest_1", 1);
//    std::cout << bl.data.size() << ' ' << bl.data_block_sizes.size() << ' ' << bl.data_block_sizes[0] << '\n';
//    double x0 = bl.data[0].x/100.;
//    double y0 = bl.data[0].y/100.;
//    std::cout << x0 << ' ' << y0 << ' ' << bl.get_name(bl.data[0].id) << ' ' << bl.data[0].energy << '\n';
//    bl.print_data("F://1.New//Corsika//output_data.txt");

    //bl.read_file_comp("W:\\DATA\\output_data_corsika.txt");
    int number_of_event = 4;

//    for(int i = number_of_event; i < 6; ++i){
    bl.read_event_comp("W:\\DATA\\output_data_corsika.txt", number_of_event);
    std::cout << "data size: " << bl.data.size() << '\n';
    //std::vector<float> core{7.8, 6.3};
    //bl.set_core(core);
    bl.clear_all();
    //bl.read_event_comp("C://Datas//output_data.txt", 1);
    //std::cout << "data size: " << bl.data.size() << '\n';

    //std::vector<float> core1{2.8, 1.3};
    //bl.set_core(core1);

    //std::vector<float> bar;
    //float num{};
    //std::ifstream input("C:\\Users\\ASUS Zephyrus\\Desktop\\ENDA\\build_noqt\\out_.txt");

    //while(input >> num){
    //    bar.push_back(num);
    //}
    //std::cout << bar.size() << '\n';

    //bl.params_size();
    //bl.conv_to_format(bar);
    //int num_of_event = 2;
    //int numberOfEvent = bl.data.size();
    //std::cout << numberOfEvent << '\n';
    //std::cout << bl.get_n_events() << '\n';
    //for (int i = bl.data.size()-10; i < bl.data.size(); ++i) std::cout << bl.data[i].id << ' ' << bl.data[i].x << ' ' << bl.data[i].y << ' ' << bl.data[i].energy << '\n';
//    bl.clear_all();
//    }
    return 0;
}
