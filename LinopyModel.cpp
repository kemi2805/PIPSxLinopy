#include "LinopyModel.hpp"


Linopy::Linopy() {

}
Linopy::Linopy(std::string filepath, int local_id, int id_size) {
    Filepath = filepath;
    id = local_id;
    nScenario = id_size;
    Set_equality_Vector();
    Set_inequality_Vector();
    Set_Matrix();
    Check_Matrix_if_empty();
}




//Code for getting the Vectors
//-------------------------------------
void Linopy::Set_equality_Vector() {
    b = GetDataFromFileEqualityConstraints<double>("b");
    c = GetDataFromFileEqualityConstraints<double>("c");
    GetDataForxvec();
    if(!id) { // Ich will das bl in Block 0 steht. Für alle anderen Blöcke ist nScenario 0.
        bl = GetDataFromFileEqualityConstraints<double>("b",true);
    }
    return;
}


void Linopy::Set_inequality_Vector() {
    auto [xl_temp,ixl_temp] = GetDataFromFileInequalityConstraints("xl");
    auto [xu_temp,ixu_temp] = GetDataFromFileInequalityConstraints("xu");
    auto [dl_temp,idl_temp] = GetDataFromFileInequalityConstraints("dl");
    auto [du_temp,idu_temp] = GetDataFromFileInequalityConstraints("du");
    if(!id) {
        auto [dll_temp,idll_temp] = GetDataFromFileInequalityConstraints("dl",true);
        auto [dlu_temp,idlu_temp] = GetDataFromFileInequalityConstraints("du",true);
        dll = std::move(dll_temp); idll = std::move(idll_temp);
        dlu = std::move(dlu_temp); idlu = std::move(idlu_temp);
    }
    xl = std::move(xl_temp); ixl = std::move(ixl_temp);
    xu = std::move(xu_temp); ixu = std::move(ixu_temp);
    dl = std::move(dl_temp); idl = std::move(idl_temp);
    du = std::move(du_temp); idu = std::move(idu_temp);
    return;    
}
//-------------------------------------

void Linopy::SetPathToFile(int id, std::string name) {
    PathToFile = Filepath;
    PathToFile += "/block";
    PathToFile += std::to_string(id);
    PathToFile += "/";
    PathToFile += name;
}

std::pair<std::unique_ptr<std::vector<double>>,std::unique_ptr<std::vector<double>>> Linopy::GetDataFromFileInequalityConstraints(std::string name, bool linking) {
    std::vector<double> Temp_Data;
    std::vector<double> Temp_Data_Active;
    std::string Fileline;
    std::ifstream File;
    if(linking)    
        SetPathToFile(nScenario+1, name);
    else
        SetPathToFile(id, name);  
    File.open(PathToFile);
    while (getline(File,Fileline)) {
        if(Fileline != "inf" && Fileline != "-inf") {
            Temp_Data.push_back(stod(Fileline));
            Temp_Data_Active.push_back(1);
        }
        else {
            Temp_Data.push_back(0);
            Temp_Data_Active.push_back(0);
        }
    }
    File.close();
    return std::make_pair(std::move(std::make_unique<std::vector<double>>(Temp_Data)),std::move(std::make_unique<std::vector<double>>(Temp_Data_Active)));
}
template<typename T>
std::unique_ptr<std::vector<T>> Linopy::GetDataFromFileEqualityConstraints(std::string name, bool linking) {
    std::vector<T> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    if(linking)    
        SetPathToFile(nScenario+1, name);
    else
        SetPathToFile(id, name);   
    File.open(PathToFile);
    while (getline(File,Fileline)) {
        Temp_Data.push_back(stod(Fileline));
    }
    if(name == "c" && Fileline.empty())
        Temp_Data.push_back(0.0);
    File.close();
    return std::move(std::make_unique<std::vector<T>>(Temp_Data));  
}

void Linopy::Set_Matrix() {
    A_row = GetDataFromRowFile("A");
    A_col = GetDataFromColFile("A");
    A_data = GetDataFromDataFile("A");

    B_row = GetDataFromRowFile("B");
    B_col = GetDataFromColFile("B");
    B_data = GetDataFromDataFile("B");

    BL_row = GetDataFromRowFile("BL");
    BL_col = GetDataFromColFile("BL");
    BL_data = GetDataFromDataFile("BL");

    C_row = GetDataFromRowFile("C");
    C_col = GetDataFromColFile("C");
    C_data = GetDataFromDataFile("C");

    D_row = GetDataFromRowFile("D");
    D_col = GetDataFromColFile("D");
    D_data = GetDataFromDataFile("D");

    DL_row = GetDataFromRowFile("DL");
    DL_col = GetDataFromColFile("DL");
    DL_data = GetDataFromDataFile("DL");

    // I decided to do it like that, so if one wants the other Format, he can just grab that;
    A_row = Transform_Matrix_Rows(std::move(A_row));
    B_row = Transform_Matrix_Rows(std::move(B_row));
    BL_row = Transform_Matrix_Rows(std::move(BL_row));
    C_row = Transform_Matrix_Rows(std::move(C_row));
    D_row = Transform_Matrix_Rows(std::move(D_row));
    DL_row = Transform_Matrix_Rows(std::move(DL_row));
    return;
}

std::unique_ptr<std::vector<long long>> Linopy::GetDataFromRowFile(std::string name, bool linking) {
    std::vector<long long> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    PathToFile += "_row";
    File.open(PathToFile);

    getline(File,Fileline);
    if(!(Fileline.size() != 0))  
        return std::move(std::make_unique<std::vector<long long>>(Temp_Data));
    Temp_Data.push_back(stoi(Fileline));
    while (getline(File,Fileline)) 
        Temp_Data.push_back(stoi(Fileline));
    File.close();
    return std::move(std::make_unique<std::vector<long long>>(Temp_Data));
}

std::unique_ptr<std::vector<long long>> Linopy::GetDataFromColFile(std::string name, bool linking) {
    std::vector<long long> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    PathToFile += "_col";
    File.open(PathToFile);
    while (getline(File,Fileline))
        Temp_Data.push_back(stoi(Fileline));
    File.close();
    return std::move(std::make_unique<std::vector<long long>>(Temp_Data));
}

std::unique_ptr<std::vector<double>> Linopy::GetDataFromDataFile(std::string name, bool linking) {
    std::vector<double> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    PathToFile += "_data";
    File.open(PathToFile);
    while (getline(File,Fileline))
        Temp_Data.push_back(stod(Fileline));
    File.close();
    return std::move(std::make_unique<std::vector<double>>(Temp_Data));
}

void Linopy::GetDataForxvec(std::string name, bool linking) { //Schreib das noch zuende. Du willst nur die Abstände wissen.
    std::vector<long long> Temp_Storage;
    long long Temp_Data;
    long long Compare_Data = 0;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    File.open(PathToFile);
    getline(File,Fileline);
    if(Fileline.empty()) {
         xvec = std::move(std::make_unique<std::vector<long long>>(Temp_Storage));
         return;
    }
    Temp_Storage.push_back(stoi(Fileline));
    Compare_Data = Temp_Storage.at(0);
    while (getline(File,Fileline)) {
        Temp_Data = stoi(Fileline);
        if ((Compare_Data + 1) == Temp_Data)
            Compare_Data++;
        else {
            Temp_Storage.push_back(Compare_Data);
            Temp_Storage.push_back(Temp_Data);
            Compare_Data = Temp_Data;
        }
    }
    File.close();
    Temp_Storage.push_back(Temp_Data);
    xvec = std::move(std::make_unique<std::vector<long long>>(Temp_Storage));
    return; 
}


void Linopy::Set_xvec_A(long long* xvec_A0, long long size) {
    std::vector<long long> Temp_Data;
    for(long long i=0; i < size ;i++) { // Setzt den Pointer von xvec_A von egal welchen Block, auf den von Block 0
        Temp_Data.push_back(*(xvec_A0+i));
    }
    xvec_A = std::move(std::make_unique<std::vector<long long>>(Temp_Data));
    return;
}





int Linopy::GetId() {
    return id;
}
void Linopy::SetId(int Id_number) {
    id = Id_number;
}

long long Linopy::Get_n() const {
    return c->size();
}

long long Linopy::Get_my() const {
    return b->size();
}

long long Linopy::Get_mz() const {
    return dl->size();
}

long long Linopy::Get_myl() const {
    return bl->size();
}

long long Linopy::Get_mzl() const {
    return dll->size();
}

long long Linopy::Get_nnz_A() const {
    return A_data->size();
}

long long Linopy::Get_nnz_B() const {
    return B_data->size();
}

long long Linopy::Get_nnz_BL() const {
    return BL_data->size();
}

long long Linopy::Get_nnz_C() const {
    return C_data->size();
}

long long Linopy::Get_nnz_D() const {
    return D_data->size();
}

long long Linopy::Get_nnz_DL() const {
    return DL_data->size();
}

long long* Linopy::Get_xvec_A() {
    return xvec->data();
}
long long Linopy::Get_xvec_A_size() {
    return xvec->size();
}



// Functions Linopy_Init calls
std::unique_ptr<std::vector<double>> Linopy::Get_equality_vector() {
    return std::move((b));
}
std::unique_ptr<std::vector<double>> Linopy::Get_equality_linking_vector() {
    return std::move(bl);
}
std::unique_ptr<std::vector<double>> Linopy::Get_objective_function() {
    return std::move(c);
}

std::unique_ptr<std::vector<double>> Linopy::Get_upper_inequality_constraint_vector() { // du
    return std::move(du);
}
std::unique_ptr<std::vector<double>> Linopy::Get_lower_inequality_constraint_vector() {
    return std::move(dl);
}
std::unique_ptr<std::vector<double>> Linopy::Get_upper_inequality_constraint_linking_vector() {
    return std::move(dlu);
}
std::unique_ptr<std::vector<double>> Linopy::Get_lower_inequality_constraint_linking_vector() {
    return std::move(dll);
}

std::unique_ptr<std::vector<double>> Linopy::Get_upper_inequality_constraint_vector_indicator() { // du
    return std::move(idu);
}
std::unique_ptr<std::vector<double>> Linopy::Get_lower_inequality_constraint_vector_indicator() {
    return std::move(idl);
}
std::unique_ptr<std::vector<double>> Linopy::Get_upper_inequality_constraint_linking_vector_indicator() {
    return std::move(idlu);
}
std::unique_ptr<std::vector<double>> Linopy::Get_lower_inequality_constraint_linking_vector_indicator() {
    return std::move(idll);
}

std::unique_ptr<std::vector<double>> Linopy::Get_upper_inequality_vector() {
    return std::move(xu);
}
std::unique_ptr<std::vector<double>> Linopy::Get_lower_inequality_vector() {
    return std::move(xl);
}
std::unique_ptr<std::vector<double>> Linopy::Get_upper_inequality_vector_indicator() {
    return std::move(ixu);
}
std::unique_ptr<std::vector<double>> Linopy::Get_lower_inequality_vector_indicator() {
    return std::move(ixl);
}

std::unique_ptr<std::vector<long long>> Linopy::Get_A_row() {
    return std::move(A_row);
}
std::unique_ptr<std::vector<long long>> Linopy::Get_A_col() {
    return std::move(A_col);
}
std::unique_ptr<std::vector<double>> Linopy::Get_A_data() {
    return std::move(A_data);
}

std::unique_ptr<std::vector<long long>> Linopy::Get_B_row() {
    return std::move(B_row);
}
std::unique_ptr<std::vector<long long>> Linopy::Get_B_col() {
    return std::move(B_col);
}
std::unique_ptr<std::vector<double>> Linopy::Get_B_data() {
    return std::move(B_data);
}

std::unique_ptr<std::vector<long long>> Linopy::Get_BL_row() {
    return std::move(BL_row);
}
std::unique_ptr<std::vector<long long>> Linopy::Get_BL_col() {
    return std::move(BL_col);
}
std::unique_ptr<std::vector<double>> Linopy::Get_BL_data() {
    return std::move(BL_data);
}

std::unique_ptr<std::vector<long long>> Linopy::Get_C_row() {
    return std::move(C_row);
}
std::unique_ptr<std::vector<long long>> Linopy::Get_C_col() {
    return std::move(C_col);
}
std::unique_ptr<std::vector<double>> Linopy::Get_C_data() {
    return std::move(C_data);
}

std::unique_ptr<std::vector<long long>> Linopy::Get_D_row() {
    return std::move(D_row);
}
std::unique_ptr<std::vector<long long>> Linopy::Get_D_col() {
    return std::move(D_col);
}
std::unique_ptr<std::vector<double>> Linopy::Get_D_data() {
    return std::move(D_data);
}

std::unique_ptr<std::vector<long long>> Linopy::Get_DL_row() {
    return std::move(DL_row);
}
std::unique_ptr<std::vector<long long>> Linopy::Get_DL_col() {
    return std::move(DL_col);
}
std::unique_ptr<std::vector<double>> Linopy::Get_DL_data() {
    return std::move(DL_data);
}

void Linopy::Transform_Matrix_Cols() {
    if (xvec_A != nullptr || id==0) {
        Transform_Matrix_Cols_A(A_col->begin(), A_col->end());
        Transform_Matrix_Cols_A(C_col->begin(), C_col->end());
    }
    if (xvec != nullptr) {
        Transform_Matrix_Cols_B(B_col->begin(), B_col->end());
        Transform_Matrix_Cols_B(BL_col->begin(), BL_col->end());
        Transform_Matrix_Cols_B(D_col->begin(), D_col->end());
        Transform_Matrix_Cols_B(DL_col->begin(), DL_col->end());
    }
    return;

}

void Linopy::Transform_Matrix_Cols_A(std::vector<long long>::iterator Col_begin, std::vector<long long>::iterator Col_end) {
    int counter = 0;
    for(auto line = Col_begin; line != Col_end; line++) {
        for(std::vector<long long>::iterator it = xvec_A->begin(); it <= xvec_A->end(); it++) {
            it++;
            if((*line <= *it) && (*line >= *(it-1))) {
                for(int i = counter; i!=-1; i--) {
                    if(!(i%2))
                        *line -= xvec_A->at(i);
                    else
                        *line += xvec_A->at(i);
                }
                *line += counter/2;
                break;
            }
            else if(it == xvec_A->end()-1) { // Weil wir in Python 0er dazufügen
                *line = 0;
                counter = 0;
                break;
            }
            counter += 2; 
        }
        counter = 0;
    }
}

void Linopy::Transform_Matrix_Cols_B(std::vector<long long>::iterator Col_begin, std::vector<long long>::iterator Col_end) {
    int counter = 0;
    for(auto line = Col_begin; line != Col_end; line++) {
        for(auto it = xvec->begin(); it != xvec->end(); it++) {
            it++;
            if((*line <= *it) && (*line >= *(it-1))) {
                for(int i = counter; i!=-1;i--) {
                    if(!(i%2))
                        *line -= xvec->at(i);
                    else
                        *line += xvec->at(i);
                }
                *line += counter/2;
                break;
            }
            else if(it == xvec->end()-1) { // Weil wir in Python 0er dazufügen
                *line = 0;
                counter = 0;
                break;
            }
            counter += 2; 
        }
        counter = 0;
    }
}

std::unique_ptr<std::vector<long long>> Linopy::Transform_Matrix_Rows(std::unique_ptr<std::vector<long long>> M_Row) {
    std::vector<long long> Temp_Data;
    long long temp = 0;
    Temp_Data.push_back(temp);
    if(M_Row->empty()) {
        Temp_Data.push_back(0);
        return std::move(std::make_unique<std::vector<long long>>(Temp_Data));
    }
    for(std::vector<long long>::iterator it = M_Row->begin() + 1; it != M_Row->end(); it++) {
        if(*(it-1) != *it) {
            temp++;
            Temp_Data.push_back(temp);
        }
        else
            temp++;
    }
    temp++;
    Temp_Data.push_back(temp);
    return std::move(std::make_unique<std::vector<long long>>(Temp_Data));
}


void Linopy::Check_Matrix_if_empty() {
    for(auto it = A_data->begin(); it != A_data->end(); it++) {
        if(*it != 0.0)
            break;
        else if(it+1 == A_data->end()) {
            A_row->clear();
            A_data->clear();
            A_col->clear();
        }
    }
        for(auto it = B_data->begin(); it != B_data->end(); it++) {
        if(*it != 0.0) 
            break;
        else if(it+1 == B_data->end()) {
            B_row->clear();
            B_data->clear();
            B_col->clear();
        }
    }
    for(auto it = BL_data->begin(); it != BL_data->end(); it++) {
        if(*it != 0.0) 
            break;
        else if(it+1 == BL_data->end()) {
            BL_row->clear();
            BL_data->clear();
            BL_col->clear();
        }
    }
    for(auto it = C_data->begin(); it != C_data->end(); it++) {
        if(*it != 0.0) {
            break;}
        else if(it+1 == C_data->end()) {
            C_row->clear();
            C_data->clear();
            C_col->clear();
        }
    }
    for(auto it = D_data->begin(); it != D_data->end(); it++) {
        if(*it != 0.0) 
            break;
        else if(it+1 == D_data->end()) {
            D_row->clear();
            D_data->clear();
            D_col->clear();
        }
    }
    for(auto it = DL_data->begin(); it != DL_data->end(); it++) {
        if(*it!=0.0) 
            break;
        else if(it+1 == DL_data->end()) {
            DL_row->clear();
            DL_data->clear();
            DL_col->clear();
        }
    }
}

void Linopy::Check_Matrix(std::unique_ptr<std::vector<long long>> Row, std::unique_ptr<std::vector<long long>> Col, std::unique_ptr<std::vector<double>> Data) {
    for(auto it = Data->begin(); it != Data->end(); it++) {
        if(*it!=0.0) 
            break;
        else if(it+1 == Data->end()) {
            
        }
    }
}


/*std::unique_ptr<std::vector<int>> Linopy::Combine_yvec(std::unique_ptr<std::vector<int>>yvec1,std::unique_ptr<std::vector<int>>yvec2){
    if(yvec1 == nullptr  && yvec2 == nullptr)
        return std::move(yvec2);
    else if(!(yvec2==nullptr) && yvec1 == nullptr)    
        return std::move(yvec2);
    else if (yvec2==nullptr && !(yvec1 == nullptr))
        return std::move(yvec1);
    else {
        std::vector<std::vector<int>> Temp_AB; 
        std::vector<int> Temp_vector;
        std::vector<std::vector<int>> sorted_Temp_AB; 
        std::vector<int> Temp_move;
        bool deleter = false;
        for(long unsigned int i = 0; i < yvec1->size(); i+=2) {
            Temp_vector.at(0) = yvec1->at(i);
            Temp_vector.at(1) = yvec1->at(i+1);
            Temp_AB.push_back(Temp_vector);
        }
        for(long unsigned int i = 0; i < yvec1->size(); i+=2) {
            Temp_vector.at(0) = yvec1->at(i);
            Temp_vector.at(1) = yvec1->at(i+1);
            Temp_AB.push_back(Temp_vector);
        }
        std::sort(Temp_AB.begin(),Temp_AB.end());
        for(long unsigned int i = 0; i < Temp_AB.size(); i++) {
            if(Temp_AB.at(i).at(1) < Temp_AB.at(i+1).at(0)) {
                sorted_Temp_AB.push_back(Temp_AB.at(i));
                deleter = false;
            }
            else if(Temp_AB.at(i).at(1) == Temp_AB.at(i+1).at(0)) {
                if(deleter) {
                    Temp_vector.at(1) = Temp_AB.at(i+1).at(1);
                    deleter = false;
                }
                else {
                    Temp_vector.at(0) = Temp_AB.at(i).at(0);
                    Temp_vector.at(1) = Temp_AB.at(i+1).at(1);
                }
                sorted_Temp_AB.push_back(Temp_vector);
            }
            else {
                if(deleter) {
                    Temp_vector.at(1) = Temp_AB.at(i+1).at(1);
                }
                else {
                    Temp_vector.at(0) = Temp_AB.at(i).at(0);
                    Temp_vector.at(1) = Temp_AB.at(i+1).at(1);
                    sorted_Temp_AB.push_back(Temp_vector);
                    deleter = true;
                }
            }
        }
        for(std::vector<int>::iterator it = sorted_Temp_AB.at(0).begin(); it != sorted_Temp_AB.back().end();it++) {
            Temp_move.push_back(*it);
        }
        return std::move(std::make_unique<std::vector<int>>(Temp_move));
    }
}*/