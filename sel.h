void createLocalA(Matrix &A,mesh m){
    float t = m.getParameter(TAU);
    A.at(0).at(0) += -t/8;  A.at(0).at(1) += t/8;
    A.at(1).at(0) += -t/8;  A.at(1).at(1) += t/8;
}

void createLocalG(Matrix &G,mesh m){
    float k = m.getParameter(K);
    float l = m.getParameter(ELEMENT_LENGTH);
    G.at(0).at(0) += k/l;  G.at(0).at(1) += -k/l;
    G.at(1).at(0) += -k/l;  G.at(1).at(1) += k/l;
}

void createLocalC(Matrix &C,mesh m){
    float lambda = m.getParameter(LAMBDA);
    float cte = 2*(lambda/6);
    C.at(0).at(0) = -cte;
    C.at(0).at(1) = cte;
    C.at(1).at(0) = -cte;
    C.at(1).at(1) = cte;
}

void createLocalH(Matrix &H,mesh m){
    float v = m.getParameter(V);
    float l = m.getParameter(ELEMENT_LENGTH);
    H.at(0).at(0) += v/l;  H.at(0).at(1) += -v/l;
    H.at(1).at(0) += -v/l;  H.at(1).at(1) += v/l;
}

void createLocalE(Matrix &E,mesh m){
    float a = m.getParameter(ALPHA);
    float k = 3 * a / 2;
    E.at(0).at(0) += -k;  E.at(0).at(1) += k;
    E.at(1).at(0) += -k;  E.at(1).at(1) += k;
}

void createLocalF(Matrix &F,mesh m){
    float delta = m.getParameter(DELTA);
    F.at(0).at(0) += -delta/2;  F.at(0).at(1) += delta/2;
    F.at(1).at(0) += -delta/2;  F.at(1).at(1) += delta/2;
}

Matrix createLocalK(int element,mesh &m){
    Matrix A,G,C,H,E,F,K;

    zeroes(A,2);
    zeroes(G,2);
    zeroes(C,2);
    zeroes(H,2);
    zeroes(E,2);
    zeroes(F,2);
    createLocalA(A,m);
    createLocalG(G,m);
    createLocalC(C,m);
    createLocalH(H,m);
    createLocalE(E,m);
    createLocalF(F,m);

    Vector row1, row2, row3, row4;


    row1.push_back(A.at(0).at(0)+G.at(0).at(0));
    row1.push_back(A.at(0).at(1)+G.at(0).at(1));
    row1.push_back(C.at(0).at(0)+ H.at(0).at(0));
    row1.push_back(C.at(0).at(1)+ H.at(0).at(1));

    row2.push_back(A.at(1).at(0)+G.at(1).at(0));
    row2.push_back(A.at(1).at(1)+G.at(1).at(1));
    row2.push_back(C.at(1).at(0)+H.at(1).at(0));
    row2.push_back(C.at(1).at(1)+H.at(1).at(1));

    row3.push_back(E.at(0).at(0));
    row3.push_back(E.at(0).at(1));
    row3.push_back(F.at(0).at(0));
    row3.push_back(F.at(0).at(1));

    row4.push_back(E.at(1).at(0));
    row4.push_back(E.at(1).at(1));
    row4.push_back(F.at(1).at(0));
    row4.push_back(F.at(1).at(1));

    K.push_back(row1);
    K.push_back(row2);
    K.push_back(row3);
    K.push_back(row4);

    return K;
}

Vector createLocalb(int element,mesh &m){
    Vector b;

    float f = m.getParameter(DELTA), l = m.getParameter(ELEMENT_LENGTH);

    b.push_back(f*l/2);
    b.push_back(f*l/2);
    b.push_back(0);
    b.push_back(0);

    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;
    int index3 = index1 + nnodes;
    int index4 = index2 + nnodes;

    K.at(index1).at(index1) += localK.at(0).at(0);
    K.at(index1).at(index2) += localK.at(0).at(1);
    K.at(index2).at(index1) += localK.at(1).at(0);
    K.at(index2).at(index2) += localK.at(1).at(1);

    K.at(index1).at(index3) += localK.at(0).at(2);
    K.at(index1).at(index4) += localK.at(0).at(3);
    K.at(index2).at(index3) += localK.at(1).at(2);
    K.at(index2).at(index4) += localK.at(1).at(3);

    K.at(index3).at(index1) += localK.at(2).at(0);
    K.at(index3).at(index2) += localK.at(2).at(1);
    K.at(index4).at(index1) += localK.at(3).at(0);
    K.at(index4).at(index2) += localK.at(3).at(1);

    K.at(index3).at(index3) += localK.at(2).at(2);
    K.at(index3).at(index4) += localK.at(2).at(3);
    K.at(index4).at(index3) += localK.at(3).at(2);
    K.at(index4).at(index4) += localK.at(3).at(3);

}

void assemblyb(element e,Vector localb,Vector &b){
    int index1 = e.getNode1() - 1;
    int index2 = e.getNode2() - 1;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
}

void ensamblaje(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs,Matrix &K,Vector &b){
    int nnodes = m.getSize(NODES);
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        element e = m.getElement(i);
        assemblyK(e,localKs.at(i),K,nnodes);
        assemblyb(e,localbs.at(i),b);
    }
}


void applyDirichlet(mesh &m,Matrix &K,Vector &b){
    for(int i=0;i<m.getSize(DIRICHLET);i++){
        condition c = m.getCondition(i,DIRICHLET);

        int index = c.getNode1()-1;

        K.erase(K.begin()+index);
        b.erase(b.begin()+index);

        for(int row=0;row<K.size();row++){
            float cell = K.at(row).at(index);
            K.at(row).erase(K.at(row).begin()+index);
            b.at(row) += -1*c.getValue()*cell;
        }
    }
}


void calculate(Matrix &K, Vector &b, Vector &T){
    cout << "Iniciando calculo de respuesta...\n";
    Matrix Kinv;
    cout << "Calculo de inversa...\n";
    inverseMatrix(K,Kinv);
    cout << "Calculo de respuesta...\n";
    productMatrixVector(Kinv,b,T);
}
