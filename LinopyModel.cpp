#include "LinopyModel.hpp"


Linopy::Linopy() {

}
Linopy::Linopy(std::string filepath, int local_id, int id_size) {
    Filepath = filepath;
    id = local_id;
    nScenario = id_size;
    std::cout << "Hier kreieren wir Linopy" << std::endl;
    Set_equality_Vector();
    std::cout << "Set_equality_vector" << std::endl;
    Set_inequality_Vector();
    std::cout << "Set_inequality_Vector" << std::endl;
    Set_Matrix();
    
    std::cout << "Set_Matrix" << std::endl;
/*  if((yvec_A->at(0) == yvec_B->at(0)) && (yvec_A->size() == yvec_B->size()))
        std::cout << "Alles läuft wie geschnitten Brot" << std::endl;
    if((yvec_C->at(0) == yvec_C->at(0)) && (yvec_D->size() == yvec_B->size()))
        std::cout << "Alles läuft wie geschnitten Brot" << std::endl;*/
}




//Code for getting the Vectors
//-------------------------------------
void Linopy::Set_equality_Vector() {
    std::cout << "In Set_equality_Vector" << std::endl;
    b = GetDataFromFileEqualityConstraints<double>("b");
    std::cout << "In Set_equality_Vector_Get_b" << std::endl;
    c = GetDataFromFileEqualityConstraints<double>("c");
    std::cout << "In Set_equality_Vector_Get_c" << std::endl;
    GetDataForxvec();
    std::cout << "Hier passiert wohl was, was nicht sollte" <<std::endl;
    if(nScenario) { // Ich will das bl in Block 0 steht. Für alle anderen Blöcke ist nScenario 0.
        bl = GetDataFromFileEqualityConstraints<double>("bl",true);
        std::cout << "In Set_equality_Vector_Get_bl" << std::endl;
    }
    std::cout << "In Set_equality_Vector_at the end" << std::endl;
    return;
}


void Linopy::Set_inequality_Vector() {
    auto [xl_temp,ixl_temp] = GetDataFromFileInequalityConstraints("xl");
    auto [xu_temp,ixu_temp] = GetDataFromFileInequalityConstraints("xu");
    auto [dl_temp,idl_temp] = GetDataFromFileInequalityConstraints("dl");
    auto [du_temp,idu_temp] = GetDataFromFileInequalityConstraints("du");
    if(nScenario) {
        auto [dll_temp,idll_temp] = GetDataFromFileInequalityConstraints("dll",true);
        auto [dlu_temp,idlu_temp] = GetDataFromFileInequalityConstraints("dlu",true);
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
        id = nScenario + 1;
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
    if(linking)    
        id = 0;
    return std::make_pair(std::move(std::make_unique<std::vector<double>>(Temp_Data)),std::move(std::make_unique<std::vector<double>>(Temp_Data_Active)));
}
template<typename T>
std::unique_ptr<std::vector<T>> Linopy::GetDataFromFileEqualityConstraints(std::string name, bool linking) {
    std::vector<T> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    if(linking)    
        id = nScenario + 1; 
    SetPathToFile(id, name);
    File.open(PathToFile);
    while (getline(File,Fileline)) {
        Temp_Data.push_back(stod(Fileline));
    }
    File.close();
    if(linking)    
        id = 0; 
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

std::unique_ptr<std::vector<int>> Linopy::GetDataFromRowFile(std::string name, bool linking) {
    std::vector<int> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    PathToFile += "_row";
    File.open(PathToFile);

    getline(File,Fileline);
    if(!(Fileline.size() != 0))  
        return std::move(std::make_unique<std::vector<int>>(Temp_Data));
    Temp_Data.push_back(stoi(Fileline));
    while (getline(File,Fileline)) 
        Temp_Data.push_back(stoi(Fileline));
    File.close();
    return std::move(std::make_unique<std::vector<int>>(Temp_Data));
}

std::unique_ptr<std::vector<int>> Linopy::GetDataFromColFile(std::string name, bool linking) {
    std::vector<int> Temp_Data;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    PathToFile += "_col";
    File.open(PathToFile);
    while (getline(File,Fileline))
        Temp_Data.push_back(stoi(Fileline));
    File.close();
    return std::move(std::make_unique<std::vector<int>>(Temp_Data));
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
    std::vector<int> Temp_Storage;
    int Temp_Data;
    int Compare_Data = 0;
    std::string Fileline;
    std::ifstream File;
    SetPathToFile(id, name);
    File.open(PathToFile);
    getline(File,Fileline);
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
    xvec = std::move(std::make_unique<std::vector<int>>(Temp_Storage));
    return; 
}


void Linopy::Set_xvec_A(int* xvec_A0, int size) {
    std::vector<int> Temp_Data;
    for(int i=0; i < size ;i++) { // Setzt den Pointer von xvec_A von egal welchen Block, auf den von Block 0
        Temp_Data.push_back(*(xvec_A0+i));
    }
    xvec_A = std::move(std::make_unique<std::vector<int>>(Temp_Data));
    return;
}





int Linopy::GetId() {
    return id;
}
void Linopy::SetId(int Id_number) {
    id = Id_number;
}

int Linopy::Get_n() const {
    return xvec->size();
}

int Linopy::Get_my() const {
    return b->size();
}

int Linopy::Get_mz() const {
    return dl->size();
}

int Linopy::Get_myl() const {
    return bl->size();
}

int Linopy::Get_mzl() const {
    return dll->size();
}

int Linopy::Get_nnz_A() const {
    return A_data->size();
}

int Linopy::Get_nnz_B() const {
    return B_data->size();
}

int Linopy::Get_nnz_BL() const {
    return BL_data->size();
}

int Linopy::Get_nnz_C() const {
    return C_data->size();
}

int Linopy::Get_nnz_D() const {
    return D_data->size();
}

int Linopy::Get_nnz_DL() const {
    return DL_data->size();
}

int* Linopy::Get_xvec_A() {
    return xvec->data();
}



// Functions Linopy_Init calls
double* Linopy::Get_equality_vector() const {
    return b->data();
}
double* Linopy::Get_equality_linking_vector() const {
    return bl->data();
}
double* Linopy::Get_objective_function() const {
    return c->data();
}

double* Linopy::Get_upper_inequality_constraint_vector() const { // du
    return du->data();
}
double* Linopy::Get_lower_inequality_constraint_vector() const {
    return dl->data();
}
double* Linopy::Get_upper_inequality_constraint_linking_vector() const {
    return dlu->data();
}
double* Linopy::Get_lower_inequality_constraint_linking_vector() const {
    return dll->data();
}

double* Linopy::Get_upper_inequality_constraint_vector_indicator() const { // du
    return idu->data();
}
double* Linopy::Get_lower_inequality_constraint_vector_indicator() const {
    return idl->data();
}
double* Linopy::Get_upper_inequality_constraint_linking_vector_indicator() const {
    return idlu->data();
}
double* Linopy::Get_lower_inequality_constraint_linking_vector_indicator() const {
    return idll->data();
}

double* Linopy::Get_upper_inequality_vector() const {
    return xu->data();
}
double* Linopy::Get_lower_inequality_vector() const {
    return xl->data();
}
double* Linopy::Get_upper_inequality_vector_indicator() const {
    return ixu->data();
}
double* Linopy::Get_lower_inequality_vector_indicator() const {
    return ixl->data();
}

int* Linopy::Get_A_row() const {
    return A_row->data();
}
int* Linopy::Get_A_col() const {
    return A_col->data();
}
double* Linopy::Get_A_data() const {
    return A_data->data();
}

int* Linopy::Get_B_row() const {
    return B_row->data();
}
int* Linopy::Get_B_col() const {
    return B_col->data();
}
double* Linopy::Get_B_data() const {
    return B_data->data();
}

int* Linopy::Get_BL_row() const {
    return BL_row->data();
}
int* Linopy::Get_BL_col() const {
    return BL_col->data();
}
double* Linopy::Get_BL_data() const {
    return BL_data->data();
}

int* Linopy::Get_C_row() const {
    return C_row->data();
}
int* Linopy::Get_C_col() const {
    return C_col->data();
}
double* Linopy::Get_C_data() const {
    return C_data->data();
}

int* Linopy::Get_D_row() const {
    return D_row->data();
}
int* Linopy::Get_D_col() const {
    return D_col->data();
}
double* Linopy::Get_D_data() const {
    return D_data->data();
}

int* Linopy::Get_DL_row() const {
    return DL_row->data();
}
int* Linopy::Get_DL_col() const {
    return DL_col->data();
}
double* Linopy::Get_DL_data() const {
    return DL_data->data();
}

void Linopy::Transform_Matrix_Cols() {
    if (xvec_A!= nullptr) {
        Transform_Matrix_Cols_A(A_col->begin(), A_col->end());
        Transform_Matrix_Cols_A(C_col->begin(), C_col->end());
    }
    if (xvec_A != nullptr) {
        Transform_Matrix_Cols_B(B_col->begin(), B_col->end());
        Transform_Matrix_Cols_B(BL_col->begin(), BL_col->end());
        Transform_Matrix_Cols_B(D_col->begin(), D_col->end());
        Transform_Matrix_Cols_B(DL_col->begin(), DL_col->end());
    }
    return;

}

void Linopy::Transform_Matrix_Cols_A(std::vector<int>::iterator Col_begin, std::vector<int>::iterator Col_end) {
    int counter = 0;
    for(auto line = Col_begin; line != Col_end; line++) {
        for(auto it = xvec_A->begin(); it != xvec_A->end(); it++) {
            it++;
            if((*line <= *it) && (*line >= *(it-1))) {
                for(int i = counter; i!=-1;i--) {
                    if(!(i%2))
                        *line -= xvec_A->at(i);
                    else
                        *line += xvec_A->at(i);
                }
                *line += counter/2;
                break;
            }
            counter += 2; 
        }
        counter = 0;
    }
}

void Linopy::Transform_Matrix_Cols_B(std::vector<int>::iterator Col_begin, std::vector<int>::iterator Col_end) {
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
            counter += 2; 
        }
        counter = 0;
    }
}

std::unique_ptr<std::vector<int>> Linopy::Transform_Matrix_Rows(std::unique_ptr<std::vector<int>> M_Row) {
    std::vector<int> Temp_Data;
    int temp = 0;
    Temp_Data.push_back(0);
    if(M_Row->empty())
        return std::move(std::make_unique<std::vector<int>>(Temp_Data));
    for(std::vector<int>::iterator it = M_Row->begin() + 1; it != M_Row->end(); it++) {
        if(*it != *(it+1)) {
            temp++;
            Temp_Data.push_back(temp);
        }
        else
            temp++;
    }
    temp++;
    Temp_Data.push_back(temp);
    return std::move(std::make_unique<std::vector<int>>(Temp_Data));
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