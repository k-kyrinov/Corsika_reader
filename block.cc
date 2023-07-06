#include "block.hh"
#include <mutex>
std::mutex m_lock;

void Block::fill_block(std::string header) {
    if (header.compare("RUNH") == 0) {
        if (run_header.size() != 0) run_header.clear();
        run_header.insert(run_header.end(), &buffer[0], &buffer[block_size]);
    }
    else if (header.compare("EVTH") == 0) {
        if (event_header.size() != 0) clear_all();
        event_header.insert(event_header.end(), &buffer[0], &buffer[block_size]);
    }
    else if (header.compare("EVTE") == 0) {
        event_end.insert(event_end.end(), &buffer[0], &buffer[block_size]);
        data_block_sizes.push_back(data.size());
    }
    else if (header.compare("RUNE") == 0) {
        run_end.insert(run_end.end(), &buffer[0], &buffer[block_size]);
    }
    else if (header.compare("LONG") == 0) {
        //..
    }
    else {
        data_block.insert(data_block.end(), &buffer[0], &buffer[block_size]);
        for (int i = 0; i <= block_size - 7; i += 7) {
            //  if((buffer[i+4]*buffer[i+4]+buffer[i+5]*buffer[i+5]) < radius_cut*radius_cut &&
       //       std::sqrt((buffer[i+1]*buffer[i+1]+buffer[i+2]*buffer[i+2]+buffer[i+3]*buffer[i+3]))*1000 /*+
       //       gParticleMass[static_cast<int>(buffer[i]/1000)]*/ > energy_cut ){

            if (buffer[i] != 0 && (buffer[i + 4] * buffer[i + 4] + buffer[i + 5] * buffer[i + 5]) < radius_cut * radius_cut &&
                std::sqrt((buffer[i + 1] * buffer[i + 1] + buffer[i + 2] * buffer[i + 2] + buffer[i + 3] * buffer[i + 3])) * 1000. > energy_cut) {
                Particle p(buffer[i], buffer[i + 1], buffer[i + 2], buffer[i + 3], buffer[i + 4], buffer[i + 5], buffer[i + 6]);
                data.push_back(p);
            }
        }
    }
}

void Block::print_header(std::string path, std::string head){
    std::ofstream output(path.c_str());
    if(head.compare("RUNH")==0){
        for(int i = 0; i < run_header.size(); ++i){
            if(i==0) output << reinterpret_cast<char*>(&run_header[0]) << '\n';
            else  output << run_header[i] << std::endl;
        }
    }
    else if (head.compare("EVTH")==0) {
        for(int i = 0; i < event_header.size(); ++i){
            if(i==0) output << reinterpret_cast<char*>(&event_header[0]) << '\n';
            else  output << event_header[i] << std::endl;
        }
    }
    else if(head.compare("EVTE")==0){
        for(int i = 0; i < event_end.size(); ++i){
            if(i==0) output << reinterpret_cast<char*>(&event_end[0]) << '\n';
            else  output << event_end[i] << std::endl;
        }
    }
    else if(head.compare("RUNE")==0){
        for(int i = 0; i < run_end.size(); ++i){
            if(i==0) output << reinterpret_cast<char*>(&run_end[0]) << '\n';
            else  output << run_end[i] << std::endl;
        }
    }
    else if(head.compare("LONG")==0){
        //..
    }
    else if(head.compare("DATA")==0){
//        for(int i = 0; i < data_block.size(); ++i){
//            output << data_block[i] << std::endl;
//        }
        for(int i = 0; i <= data_block.size()-7; i+=7){
            output << int(data_block[i]/1000) << ' ' << int(std::sqrt((data_block[i+1]*data_block[i+1]+data_block[i+2]*data_block[i+2]+data_block[i+3]*data_block[i+3]))*1000)
            << ' '<< int(std::sqrt(data_block[i+4]*data_block[i+4]+data_block[i+5]*data_block[i+5])) << ' ' << data_block[i+4]  << '\n';
        }
    }
}

// RUN HEADER
float Block::get_energy_min(){
    return run_header[16];
}
float Block::get_energy_max(){
    return run_header[17];
}
float Block::get_nevents(){
    return run_header[92];
}
float Block::rh_get_element(int index){
    return run_header[index];
}

// EVENT HEADER
float Block::get_primary_id(){
    return event_header[2];
}
float Block::get_total_energy(){
    return event_header[3];
}
float Block::get_zenith(){
    return event_header[10];
}
float Block::get_azimuth(){
    return event_header[11];
}
float Block::eh_get_element(int index){
    return event_header[index];
}

// EVENT END
float Block::get_photons(){
    return event_end[2];
}
float Block::get_electrons(){
    return event_end[3];
}
float Block::get_hadrons(){
    return event_end[4];
}
float Block::get_muons(){
    return event_end[5];
}
float Block::get_particles(){
    return event_end[6];
}
float Block::ee_get_element(int index){
    return event_end[index];
}

// RUN END
float Block::get_run_number(){
    return run_end[1];
}
float Block::get_n_events(){
    return run_end[2];
}
float Block::re_get_element(int index){
    return run_end[index];
}

// DATA BLOCK
void Block::print_data(std::string path){
    std::ofstream output(path.c_str());
    for(int i = 0; i < data.size(); ++i){
        output << data[i].id << ' ' << data[i].energy << ' '<< data[i].x << ' '<< data[i].y << ' '<< data[i].t << '\n';
    }
}

void Block::read_file(std::string filename) {
    std::ifstream file;
    //std::ofstream output("F://1.New//Corsika//test_out.txt");
    float sizeBuff;
    file.open(filename.c_str(), std::fstream::binary | std::fstream::out); //p500TeVtest_1

    int j = 0;
    if (file) {
        file.seekg(0);
        file.read((char*)&sizeBuff, 4); // junk
        while (!file.eof()) {
            file.read((char*)&buffer, sizeof(buffer));
            for (int i = 0; i < 273; ++i) {
                if (i == 0) {
                    std::string oString;
                    oString = reinterpret_cast<char*>(&buffer[0]);
                    fill_block(oString);
                }
            }

            j += 1;
            if (j % 21 == 0) {
                for (int i = 0; i < 2; ++i) {
                    file.read((char*)&sizeBuff, 4);
                }
            }
        }
    }
}

void Block::read_event(std::string filename, int event) {
    std::ifstream file;
    //std::ofstream output("F://1.New//Corsika//test_out.txt");
    float sizeBuff;
    file.open(filename.c_str(), std::fstream::binary | std::fstream::out);  /*"F://1.New//Corsika//p500TeVtest_1"*/

    int j = 0;
    //int num_of_event{};
    if (file) {
        file.seekg(0);
        file.read((char*)&sizeBuff, 4);
        while (!file.eof()) {
            file.read((char*)&buffer, sizeof(buffer));
            for (int i = 0; i < 273; ++i) {
                if (i == 0) {
                    std::string oString;
                    oString = reinterpret_cast<char*>(&buffer[0]);
                    if (oString.compare("EVTH") == 0) num_of_event += 1;
                    // if(oString.compare("RUNH")==0) num_of_event+=1;
                    if (event == num_of_event) fill_block(oString);
                    if (num_of_event > event || oString.compare("EVTE") == 0) break;
                }
            }

            j += 1;
            if (j % 21 == 0) {
                for (int i = 0; i < 2; ++i) {
                    file.read((char*)&sizeBuff, 4);
                }
            }
        }
    }
}

void Block::read_file_comp(std::string filename) {
    std::ifstream file(filename);
    std::string s;
    while (std::getline(file, s) )
    {
        std::stringstream lineStream(s);
        std::string cell;
        std::vector<std::string> parsedRow;
        while (std::getline(lineStream, cell, ' '))
        {
            parsedRow.push_back(cell);
        }
        if (parsedRow[0] == "event") total_nevents += 1;
    }
}

void Block::read_event_comp(std::string filename, int event) {
    std::ifstream file(filename);
    std::string s;
    std::vector<float> v;
    int num_ev{};
    int count = 0;
    while (std::getline(file, s))
    {
        count += 1;
        std::stringstream lineStream(s);
        std::string cell;
        std::vector<std::string> parsedRow;
        while (std::getline(lineStream, cell, ' '))
        {
            parsedRow.push_back(cell);
        }
        if (parsedRow[0] == "event") {
            ++num_ev;
            if (num_ev > event) break;
            //std::cout << parsedRow[0] << ' ' << parsedRow[1] << ' ' << parsedRow[2] << ' ' << parsedRow[3] << ' ' << parsedRow[4] << ' ' << parsedRow[5] << '\n';
            //if (params.size() != 0) params.clear();
            if (std::atoi(parsedRow[1].c_str()) == event && num_ev == event ) {
                params.insert(params.end(), { /*std::atof(parsedRow[2].c_str()),*/ std::atof(parsedRow[3].c_str()), std::atof(parsedRow[4].c_str()),
                              std::atof(parsedRow[5].c_str()), std::atof(parsedRow[6].c_str()) });
            }
        }
        else if (parsedRow[0] != "event" && num_ev == event) {
            std::transform(parsedRow.begin(), parsedRow.end(), std::back_inserter(v),
                [&](std::string s) {
                    return stoi(s);
                });
            Particle p(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7]);   //   id, energy, px, py, pz, x, y, t
            data.push_back(p);
            v.clear();
        }
    }
}

void Block::set_core(std::vector<float>& core) {
    params.insert(params.end(), { core[0], core[1]});
}

void Block::conv_to_format(std::vector<float> bar, std::string filename, int *Ns, int *N_th) {
    std::ofstream output;
    output.open(filename.c_str(), std::ios_base::app);
    std::vector<std::vector<double>> ed_det(number_of_det, std::vector<double>());
    std::vector<std::vector<double>> times_det(number_of_det, std::vector<double>());
    std::vector<std::vector<int>> num_n_det(number_of_det, std::vector<int>());
    std::vector<std::vector<int>> num_mu_det(number_of_det, std::vector<int>());
    std::string s;
    int part_id{}; // Corsika info
    double x_core{}, y_core{}, E_0{}, theta{}, phi{}; // Corsika info

    int det_Num{}, trigged_det{}, M1{}, M2{}, M{}; // N_thr - номер броска n-го ливня
    double temp{}, temp_n{}, temp_mu{}, Sum_ed{}; //, mip = 450; // keV ??
    int Num_n{}, Num_mu{};
    for (int k = 0; k < bar.size(); k += 5) {
        if (bar[k] == -1) {
            k += 1;
            for (int i = 0; i < number_of_det; ++i) {
                temp = 0;
                temp_n = 0;
                temp_mu = 0;
                for (int j = 0; j < ed_det[i].size(); ++j) {
                    temp += ed_det[i][j];
                    temp_n += num_n_det[i][j];
                    temp_mu += num_mu_det[i][j];
                    Num_n += num_n_det[i][j];
                    Num_mu += num_mu_det[i][j];
                }
                ed_det[i].clear();
                num_n_det[i].clear();
                num_mu_det[i].clear();
                num_n_det[i].push_back(temp_n);
                num_mu_det[i].push_back(temp_mu);
                ed_det[i].push_back(temp / mip);
                Sum_ed += temp / mip/ rp;
            }

            for (int i = 0; i < number_of_det; ++i) {
                if (times_det[i].size() > 0) {
                    sort(times_det[i].begin(), times_det[i].end());
                    temp = times_det[i][0];
                    times_det[i].clear();
                    times_det[i].push_back(temp);
                }
                else {
                    times_det[i].push_back(-1);
                }
            }

            for (int i = 0; i < number_of_det; ++i) {
                ed_det[i][0] /= rp;
                if (ed_det[i][0] > treshold) trigged_det += 1;
            }
            if (trigged_det >= krat) ++M1;
            if (Sum_ed >= 300) ++M2;

            if (M1 >= 1) M = 1;
            if (M2 >= 1) M += 2;
            if (Num_n >= 10) M += 4;
            if (M1 > 1 || M2 > 1) M += 8;

            x_core = params[4 + (*Ns - 1) * 6];
            y_core = params[5 + (*Ns - 1) * 6];
            part_id = params[(*Ns - 1) * 6];
            E_0 = params[1 + (*Ns - 1) * 6];
            theta = params[2 + (*Ns - 1) * 6];
            phi = params[3 + (*Ns - 1) * 6];
            
            // x, y, Event number, Number throws, Master, fold, Sum n, Sum m, [16*4] ed[i], num_n[i], num_m[i], times[i], 
            // part. id, E0, theta, phi

            if (M >= 1) {
                output << x_core << ' ' << y_core << ' ' << *Ns << ' ' << *N_th << ' ' << M << ' ' << trigged_det << ' ' << Num_n << ' ' << Num_mu << ' ';
                for (int i = 0; i < number_of_det; ++i) output << ed_det[i][0] << ' ' << num_n_det[i][0] << ' ' << num_mu_det[i][0] << ' ' << times_det[i][0] << ' ';
                output << part_id << ' ' << E_0 << ' ' << theta << ' ' << phi << '\n';
            }

            temp = temp_n = temp_mu = Num_n = Num_mu = trigged_det = Sum_ed = M1 = M2 = M = part_id = x_core = y_core = E_0 = theta = phi = 0;
            for (int i = 0; i < number_of_det; ++i) {
                ed_det[i].clear();
                num_n_det[i].clear();
                num_mu_det[i].clear();
                times_det[i].clear();
            }
        }
        else {
            det_Num = static_cast<int>(bar[k]); // endaHit->GetChamberNb()
            ed_det[det_Num].push_back(bar[k + 1]); // endaHit->GetEdep()
            num_n_det[det_Num].push_back(bar[k + 2]); // endaHit->GetNumNeutrons()
            num_mu_det[det_Num].push_back(bar[k + 3]); // endaHit->GetNumMuons()
            times_det[det_Num].push_back(bar[k + 4]); // endaHit->GetGlobalTime()
        }
    }
}

void Block::clear_all() {
    data.clear();
    //run_header.clear();
    event_header.clear();
    event_end.clear();
    run_end.clear();
    data_block.clear();
}