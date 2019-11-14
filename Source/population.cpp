/* Declaration of Population class */

# ifndef NSGAII_SOURCE_POPULATION_CPP
# define NSGAII_SOURCE_POPULATION_CPP

# include <algorithm>
# include <cstdlib>
# include <string>
# include <list>
# include <iostream>

# include "../Header/population.hpp"
# include "./rand.cpp"
# include "./log.cpp"

// Individual

void Individual::_decode(std::vector<int> &nbits, std::vector<double> &min_binvar, std::vector<double> &max_binvar){
    int nbin = nbits.size();
    if (nbin != 0){
        int sum, tempBinary;
        for (int i = 0; i < nbin; i++){
            sum = 0;
            tempBinary = 0;
            // gene最高位为1
            if (Individual::gene[i][0] == 1){
                sum += pow(2, nbits[i] - 1);
                tempBinary = 1;
            }
            // 将格雷码转换为二进制码
            for (int j = 0; j < nbits[i]; j++){
                // gene[j][k] 与 tempBinary相异（一个为0，一个为1），此时二进制该位为1
                if ((Individual::gene[i][j] + tempBinary) == 1){
                    sum += pow(2, nbits[i] - 1 - j);
                    tempBinary = 1;
                }else{
                    tempBinary = 0;
                }
            }
            Individual::xbin.push_back(min_binvar[i] + sum * (max_binvar[i] - min_binvar[i]) / (pow(2, nbits[i]) - 1));
        } // for (int i = 0; i < nbin; i++)
    } // if (nbin != 0)
} // void Individual::_decode

Individual::Individual(params p){
    Individual::power_dist = 0;
    Individual::crowd_dist = 0;
    Individual::cnt = 0;
    // 生成位于[0,1)间的随机数
    std::uniform_real_distribution<double> dis(0, 1);
    if (p.nreal != 0){
        for (int i = 0; i < p.nreal; i++){
            double num = dis(gen) * (p.max_realvar[i] - p.min_realvar[i]) + p.min_realvar[i];
            Individual::xreal.push_back(num);
        }

        // 编码得到gene

    } // if (p.nreal != 0)
    if (p.nbin != 0){
        // 初始化gene
        for (int i = 0; i < p.nbin; i++){
            std::vector<int> generow;
            for (int j = 0; j < p.nbits[i]; j++){
                double num = dis(gen);
                if (num > 0.5){
                    generow.push_back(1);
                }else{
                    generow.push_back(0);
                }
            }
            Individual::gene.push_back(generow);
        }
        // for (int i = 0; i < gene.size(); i++){
        //     for (int j = 0; j < gene[i].size(); j++){
        //         cout << gene[i][j] << ((j == gene[i].size() - 1)?"\n":", ");
        //     }
        // }
        
        // 译码gene
        Individual::_decode(p.nbits, p.min_binvar, p.max_binvar);
    } // if (p.nbin != 0)
} // Individual::Individual

// end Individual

// Population

void Population::_computeViolation(){
    int ncon = Population::ind[0]->constr.size();

    if (ncon == 0){
        return;
    }else{
        int popsize = Population::ind.size();
        for (int i = 0; i < popsize; i++){
            int violation = 0;
            for (int j = 0; j < ncon; j++){
                if (Population::ind[i]->constr[j] < 0.0){
                    violation += Population::ind[i]->constr[j];
                }
            } // for (int j = 0; j < ncon; j++)
            Population::ind[i]->constr_violation = violation;
        } // for (int i = 0; i < popsize; i++)
    } // if (ncon == 0)
} // void Population::_computeViolation(int ncon)


// Return values present dominates relationship
// 1 for a dominates b
// -1 for b dominates a
// 0 for both a and b are non_dominated 
int Check_dominance(Individual *a, Individual *b){
    if (a->constr_violation < 0 && b->constr_violation < 0){
        if (a->constr_violation > b->constr_violation){
            return 1;
        }else if(a->constr_violation < b->constr_violation){
            return 0;
        }else{
            return -1;
        } // if (a.constr_violation > b.constr_violation)
    }else{
        if (a->constr_violation < 0 && b->constr_violation == 0){
            return -1;
        }else if (a->constr_violation == 0 && b->constr_violation < 0){
            return 1;
        }else{
            int nobj = a->obj.size();
            int flag1 = 0, flag2 = 0;
            for (int i = 0; i < nobj; i++){
                if (a->obj[i] < b->obj[i]){
                    flag1 = 1;
                }else if (a->obj[i] > b->obj[i]){
                    flag2 = 1;
                } // if (a.obj[i] < b.obj[i])
            } // for (int i = 0; i < nobj; i++)
            
            if (flag1 == 1 && flag2 == 0){
                return 1;
            }else if (flag1 == 0 && flag2 == 1){
                return -1;
            }else{
                return 0;
            }
        }
    } // if (a.constr_violation < 0 && b.constr_violation < 0)
} // int Check_dominance(Individual a, Individual b)

void Population::_assign_rank_and_crowding_distance(){
    // 现在进行排序的个体的级别
    int rank = 1;
    // 准备被排序的序号
    std::list<int> cur;
    // 暂未被排序的个体的序号
    std::list<int> origin;
    int popsize = Population::ind.size();
    for (int i = 0; i < popsize; i++){
        origin.push_back(i);
    }
    // 用于定位序号的迭代器
    auto temp1 = origin.begin(), temp2 = cur.begin();

    while (!origin.empty()){
        // 未排序的元素仅剩一个的情况
        if (origin.size() == 1){
            Population::ind[origin.front()]->rank = rank;
            Population::ind[origin.front()]->crowd_dist = INF;
            break;
        }

        cur.push_front(*temp1);
        temp2 = cur.begin();
        origin.erase(temp1++);
        int front_size = 1;
        while (temp1 != origin.end()){
            temp2 = cur.begin();
            int flag = -1;

            while (temp2 != cur.end()){
                flag = Check_dominance(Population::ind[*temp1], Population::ind[*temp2]);

                if (flag == 1){
                    origin.push_front(*temp2);
                    cur.erase(temp2++);
                    front_size--;
                }else if (flag == 0){
                    temp2++;
                }else{
                    temp2 = origin.end();
                } // if (flag == 1)
            } // while (temp2 != cur.end())

            if (flag == 0 || flag == 1){
                cur.push_front(*temp1);
                front_size++;
                origin.erase(temp1++);
            }
        } // while (temp1 != origin.end())
        
        temp2 = cur.begin();
        while (temp2 != cur.end()){
            Population::ind[*temp2]->rank = rank;
            temp2++;
        } // while (temp2 != cur.end())
        
        // 原算法中assign_crowding_distance_list部分
        if (front_size < 3){
            temp2 = cur.begin();
            Population::ind[*temp2]->crowd_dist = INF;
            if (front_size == 2){
                Population::ind[*++temp2]->crowd_dist = INF;
            }

            temp2 = cur.end();
        }else{
            std::vector<int> dist;
            temp2 = cur.begin();
            for (; temp2 != cur.end(); temp2++){
                dist.push_back(*temp2);
                Population::ind[*temp2]->crowd_dist = 0.0;
            }
            // 目标函数的个数
            int nobj = Population::ind[0]->obj.size();
            // 原算法中assign_crowding_distance部分
            // 对每个目标函数进行快排
            auto it = dist.begin();
            for (int i = 0; i < nobj; i++){
                // 将dist中的序号按个体特征值从小到大排序
                std::sort(dist.begin(),dist.end(), [this, i](int a, int b)->bool{return this->ind[a]->obj[i]<this->ind[b]->obj[i];});
                
                it = dist.begin();
                Population::ind[*it]->crowd_dist = INF;
                it++;
                while (it != dist.end()){
                    // 如果拥挤距离为最大则无需更新
                    if (Population::ind[*it]->crowd_dist != INF){
                        double head = Population::ind[*dist.begin()]->obj[i], tail = Population::ind[*(dist.end()-1)]->obj[i];
                        if (head == tail){
                            Population::ind[*it] ->crowd_dist += 0.0;
                        }else{
                            Population::ind[*it]->crowd_dist += (Population::ind[*it + 1]->obj[i] + Population::ind[*it + 1]->obj[i])/(tail - head);
                        } // if (head == tail)

                        it++;
                    }else{
                        it = dist.end();
                    } // if (Population::ind[*it].crowd_dist != INF)
                } // while (it != dist.end())
            } // for (int i = 0; i < nobj; i++)
            
            it = dist.begin();
            while (it != dist.end()){
                if (Population::ind[*it]->crowd_dist != INF){
                    Population::ind[*it]->crowd_dist /= nobj;
                }

                it++;
            } // while (it != dist.end())
            // end assign_crowding_distance
        }// if(front_size < 3)
        // end assign_croding_distance_list
        cur.clear();
        rank++;
    } // while (!origin.empty())
} // Population::_assign_rank_and_crowding_distance

// 返回两个个体中处于支配地位的一方，
// 如果两个个体间不存在支配关系则随机返回一个个体
Individual* FetchInd(Individual *a, Individual *b){
    int flag;
    flag = Check_dominance(a, b);

    if (flag == 1){
        return a;
    }else if (flag == -1){
        return b;
    }else{
        if (a->crowd_dist > b->crowd_dist){
            return a;
        }else if (a->crowd_dist < b->crowd_dist){
            return b;
        }else{
            std::uniform_real_distribution<> dis(0, 1);
            if (dis(gen) > 0.5){
                return b;
            }else{
                return a;
            }
        } // if (a->crowd_dist > b->crowd_dist)
    } // if (flag == 1)
} // Individual* FetchInd(Individual *a, Individual *b)

// 接受5个参数, 两个父个体, 两个子个体
// flag =  1 意味着nreal不为0
// flag = -1 意味着nbin不为0
void Population::_crossover(Individual *p1, Individual *p2, Individual *c1, Individual *c2, int flag){
    // 生成[0.0, 1.0)间的随机数
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    double rand = dis(gen);

    if (flag == 1){
        // nreal != 0
        int nreal = this->ind[0]->xreal.size();
        double cc1, cc2;
        double alpha, beta, betaq;
        if (rand <= this->pcross_real){
            for (int i = 0; i < nreal; i++){
                rand = dis(gen);
                if (rand <= 0.5){
                    if (fabs(p1->xreal[i] - p2->xreal[i]) > EPS){
                        double y1 = 0, y2 = 0, yl = 0, yu = 0;
                        if (p1->xreal[i] < p2->xreal[i]){
                            y1 = p1->xreal[i];
                            y2 = p2->xreal[i];
                        }else{
                            y1 = p2->xreal[i];
                            y1 = p1->xreal[i];
                        }
                        yl = this->min_realvar->at(i);
                        yu = this->max_realvar->at(i);
                        rand = dis(gen);
                        beta = 1.0 + (2.0 * (y1 - yl) / (y2 - y1));
                        alpha = 2.0 - pow(beta, -(this->eta_c + 1.0));
                        if (rand <= (1.0 / alpha)){
                            betaq = pow((rand*alpha), (1.0/(this->eta_c + 1.0)));
                        }else{
                            betaq = pow((1.0 /(2.0 - rand * alpha)), (1.0 / (this->eta_c + 1.0)));
                        } // if (rand <= (1.0 / alpha))
                        cc1 = 0.5 * ((y1 + y2) - betaq * (y2 - y1));
                        beta = 1.0 + (2.0 * (yu - y2) / (y2 - y1));
                        alpha = 2.0 - pow(beta, -(this->eta_c + 1.0));
                        if (rand <= (1.0 / alpha)){
                            betaq = pow ((rand * alpha), (1.0 / (this->eta_c + 1.0)));
                        }else{
                            betaq = pow((1.0 / (2.0 - rand * alpha)), (1.0 / (this->eta_c + 1.0)));
                        } // if (rand <= (1.0 / alpha))
                        cc2 = 0.5 * ((y1 + y2) + betaq * (y2 -y1));
                        cc1 = (cc1 < yl)?yl:cc1;
                        cc1 = (cc1 > yu)?yu:cc1;
                        cc2 = (cc2 < yl)?yl:cc2;
                        cc2 = (cc2 > yu)?yu:cc2;
                        
                        if (dis(gen) <= 0.5){
                            c1->xreal[i] = cc2;
                            c2->xreal[i] = cc1;
                        }else{
                            c1->xreal[i] = cc1;
                            c2->xreal[i] = cc2;
                        } // if (dis(gen) <= 0.5)
                    } // if (fabs(p1->xreal[i] - p2->xreal[i]) > EPS)
                }else{
                    c1->xreal[i] = p1->xreal[i];
                    c2->xreal[i] = p2->xreal[i];
                } // if (rand <= 0.5)
            } // for (int i = 0; i < nreal; i++)
        }else{
            // 不进行交叉
            for (int i = 0; i < nreal; i++){
                c1->xreal[i] = p1->xreal[i];
                c2->xreal[i] = p2->xreal[i];
            } // for (int i = 0; i < nreal; i++)
        } // if (rand <= this->pcross_real)
    }else if(flag == -1){
        // nbin != 0
        int nbin = this->nbits->size();
        int site1, site2;
        for (int i = 0; i < nbin; i++){
            rand = dis(gen);
            if (rand <= this->pcross_bin){
                site1 = dis(gen) * (this->nbits->at(i) - 1);
                site2 = dis(gen) * (this->nbits->at(i) - 1);
                if (site1 > site2){
                    int temp = site1;
                    site1 = site2;
                    site2 = site1;
                }

                for (int j = 0; j < site1; j++){
                    c1->gene[i][j] = p1->gene[i][j];
                    c2->gene[i][j] = p2->gene[i][j];
                }
                for (int j = site1; j < site2; j++){
                    c1->gene[i][j] = p2->gene[i][j];
                    c2->gene[i][j] = p1->gene[i][j];
                }
                for (int j = 0; j < this->nbits->at(i); j++){
                    c1->gene[i][j] = p1->gene[i][j];
                    c2->gene[i][j] = p2->gene[i][j];
                }
            }else{
                for (int j = 0; j < this->nbits->at(i); j++){
                    c1->gene[i][j] = p1->gene[i][j];
                    c2->gene[i][j] = p2->gene[i][j];
                }
            } // if (rand <= this->pcross_bin)
        } // for (int i = 0; i < nbin; i++)
    } // if (flag == 1)
} // void Population::_crossover(Individual *p1, Individual *p2, Individual *c1, Individual *c2)

void Population::_selection(std::vector<Individual*> &childInd){
    int popsize = this->ind.size();
    // 生成[0, popsize-1]区间内的随机整数
    std::uniform_int_distribution<> dis(0, popsize - 1);

    // 初始化索引数组
    std::vector<int> index1, index2;
    for (int i = 0; i < popsize; i++){
        index1.push_back(i);
        index2.push_back(i);
    }
    int rand = 0, temp = 0;
    for (int i = 0; i < popsize; i++){
        rand = dis(gen);
        temp = index1[rand];
        index1[rand] = index1[i];
        index1[i] = temp;

        rand = dis(gen);
        temp = index2[rand];
        index2[rand] = index2[i];
        index2[i] = temp;
    }

    Individual *p1, *p2;
    for (int i = 0; i < popsize; i+=4){
        p1 = FetchInd(this->ind[index1[i]], this->ind[index1[i+1]]);
        p2 = FetchInd(this->ind[index1[i+2]], this->ind[index1[i+3]]);
    }
} // void Population::_selection()

Population::Population(params p){
    this->currentGen = 1;
    this->pcross_real = p.pcross_real;
    this->pcross_bin = p.pcross_bin;
    this->min_realvar = NULL;
    this->min_binvar = NULL;
    this->max_realvar = NULL;
    this->max_binvar = NULL;
    this->nbits = NULL;
    this->eta_c = 0;

    if (p.nreal != 0){
        this->min_realvar = new std::vector<double>(p.min_realvar);
        this->max_realvar = new std::vector<double>(p.max_realvar);
        this->eta_c = p.eta_c;
    }

    if (p.nbin != 0){
        this->nbits = new std::vector<int>(p.nbits);
        this->min_binvar = new std::vector<double>(p.min_binvar);
        this->max_binvar = new std::vector<double>(p.max_binvar);
    }

    for (int i = 0; i < p.popsize; i++){
        Individual *in = new Individual(p);
        Population::ind.push_back(in);
    }
} // Population::Population

void Population::Evolution(){
    Population::_assign_rank_and_crowding_distance();
    ReportPop(InitialPopPath);

    std::vector<Individual*> childInd;
    for (int i = 0; i < this->ind.size(); i++){
        Individual *child = new Individual(*(this->ind[i]));
    }
} // void Population::Evolution()

void Population::ReportPop(std::string path){
    ofstream file;
    file.open(LogPath, ios::app);
    if (file.is_open()){
        for (int i; i < this->ind.size(); i++){
            // 写入目标值
            int nobj = this->ind[i]->obj.size();
            file << '"';
            for (int j = 0; j < nobj; j++){
                file << this->ind[i]->obj[j] << ((j == nobj - 1)?'"':',');
            } // for (int j = 0; j < nobj; j++)
            file << ',';

            int ncon = this->ind[i]->constr.size();
            if (ncon != 0){
                file << '"';
                for (int j = 0; j < ncon; j++){
                    file << this->ind[i]->constr[j] << ((j == ncon - 1)?'"':',');
                } // for (int j = 0; j < ncon; j++)
                file << ',';
            } // if (ncon != 0)

            int nreal = this->ind[i]->xreal.size();
            if (nreal != 0){
                file << '"';
                for (int j = 0; j < nreal; j++){
                    file << this->ind[i]->xreal[j] << ((j == nreal - 1)?'"':',');
                } // for (int j = 0; j < nreal; j++)
                file << ',';
            } // if (nreal != 0)

            int nbin = this->ind[i]->gene.size();
            if (nbin != 0){
                file << '"';
                for (int j = 0; j < nbin; j++){
                    file << '"';
                    for (int k = 0; k < this->ind[i]->gene[j].size(); k++){
                        file << this->ind[i]->gene[j][i] << ((k == this->ind[i]->gene[j].size() - 1)?'"':',');
                    } // for (int k = 0; k < this->ind[i]->gene[j].size(); k++)
                    file << ',';
                }
            } // if (nbin != 0)
            file << this->ind[i]->constr_violation << ',';
            file << this->ind[i]->rank << ',';
            file << this->ind[i]->crowd_dist << endl;
        } // for (int i; i < this->ind.size(); i++)
    }else{
        RaiseError("Log file open fail, please check log path");
    } // if (file.is_open())

    file.close();
} // void Population::ReportPop(std::string path)

// end Population
# endif // NSGAII_SOURCE_POPULATION_CPP