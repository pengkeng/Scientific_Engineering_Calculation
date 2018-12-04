import Jama.Matrix;

public class Lab_2_Main {

    public void start() {

        double[][] a = new double[][]{{2, 4, 0}, {3, -1, 1}, {-2, -2, 0}};
        double[] b = new double[]{5, 9, 3};
        double[][] a1 = new double[][]{{10, -7, 0, 1}, {-3, 2.099999, 6, 2}, {-5, -1, 5, -1}, {2, 1, 0, 2}};
        double[] b1 = new double[]{8, 5.9000001, 5, 1};
        double[][] a2 = new double[][]{{2, 2, 3}, {4, 7, 7}, {-2, 4, 5}};
        double[] b2 = new double[]{3, 1, -7};


//        Gauss(new double[][]{{10, 3, 1}, {2, -10, 3}, {1, 3, 10}}, new double[]{14, -5, 15});
//        Doolittle(new double[][]{{10, 3, 1}, {2, -10, 3}, {1, 3, 10}}, new double[]{14, -5, 15});
//        jacobi(new double[][]{{10, 3, 1}, {2, -10, 3}, {1, 3, 10}}, new double[]{14, -5, 15});
//        GaussSeidel(new double[][]{{10, 3, 1}, {2, -10, 3}, {1, 3, 10}}, new double[]{14, -5, 15});
//        SuccessiveOverRelaxation(new double[][]{{10, 3, 1}, {2, -10, 3}, {1, 3, 10}}, new double[]{14, -5, 15}, 1.1);
//
//
        Gauss(getMatrixA(4), getMatrixB(4));
        Doolittle(getMatrixA(4), getMatrixB(4));
        jacobi(getMatrixA(4), getMatrixB(4));
        GaussSeidel(getMatrixA(4), getMatrixB(4));
        SuccessiveOverRelaxation(getMatrixA(4), getMatrixB(4), 1.1);


//        Gauss(new double[][]{{10, -7, 0, 1}, {-3, 2.099999, 6, 2}, {-5, -1, 5, -1}, {2, 1, 0, 2}}, new double[]{8, 5.9000001, 5, 1});
//        Doolittle(new double[][]{{10, -7, 0, 1}, {-3, 2.099999, 6, 2}, {-5, -1, 5, -1}, {2, 1, 0, 2}}, new double[]{8, 5.9000001, 5, 1});
//        jacobi(new double[][]{{10, -7, 0, 1}, {-3, 2.099999, 6, 2}, {-5, -1, 5, -1}, {2, 1, 0, 2}}, new double[]{8, 5.9000001, 5, 1});
//        GaussSeidel(new double[][]{{10, -7, 0, 1}, {-3, 2.099999, 6, 2}, {-5, -1, 5, -1}, {2, 1, 0, 2}}, new double[]{8, 5.9000001, 5, 1});
//        SuccessiveOverRelaxation(new double[][]{{10, -7, 0, 1}, {-3, 2.099999, 6, 2}, {-5, -1, 5, -1}, {2, 1, 0, 2}}, new double[]{8, 5.9000001, 5, 1}, 0.9);
//
//        ConjugateGradient(getMatrixA(4), getMatrixB(4));


    }


    private double[][] getMatrixA(int n) {
        double[][] A = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = (double) 1 / (i + 1 + j);
            }
        }
        return A;
    }

    private double[] getMatrixB(int n) {
        double[] B = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                B[i] += (double) 1 / (i + 1 + j);
            }
        }
        return B;
    }

    private double[][] cholesky(double[][] x) {
        int n = x.length;
        double[][] L1 = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                L1[i][j] = 0;
            }
        }

        L1[0][0] = Math.sqrt(x[0][0]);
        //计算第一列
        for (int i = 1; i < n; i++) {
            L1[i][0] = x[i][0] / L1[0][0];
        }

        for (int k = 1; k < n; k++) {
            double count = 0;
            for (int i = 0; i < k; i++) {
                count += Math.pow(L1[k][i], 2);
            }
            L1[k][k] = Math.sqrt(x[k][k] - count);
            for (int i = k + 1; i < n; i++) {
                double count2 = 0;
                for (int j = 0; j < k; j++) {
                    count2 += (L1[i][j] * L1[k][j]);
                }
                L1[i][k] = (x[i][k] - count2) / L1[k][k];
            }
        }
        return L1;
    }

    // 列主元高斯消元法
    private void Gauss(double a[][], double b[]) {
        int n = b.length;
        double[] x = new double[n];
        for (int k = 0; k < n - 1; k++) {
            double max = Math.abs(a[k][k]);
            int maxrow = k;
            for (int i = k; i < n; i++) {
                if (Math.abs(a[i][k]) > max) {
                    max = Math.abs(a[i][k]);
                    maxrow = i;
                }
            }
            if (max == 0) {
                return;
            }
            if (maxrow != k) {
                double temp[] = a[k];
                a[k] = a[maxrow];
                a[maxrow] = temp;
                double temp2 = b[k];
                b[k] = b[maxrow];
                b[maxrow] = temp2;
            }

            double[][] m = new double[n][n];
            for (int i = k + 1; i < n; i++) {
                m[i][k] = a[i][k] / a[k][k];
                for (int j = k + 1; j < n; j++) {
                    a[i][j] = a[i][j] - m[i][k] * a[k][j];
                }
                b[i] = b[i] - m[i][k] * b[k];
            }
        }
        x[n - 1] = b[n - 1] / a[n - 1][n - 1];
        for (int i = n - 2; i >= 0; i--) {
            double count = 0;
            for (int j = i + 1; j < n; j++) {
                count += a[i][j] * x[j];
            }
            x[i] = (b[i] - count) / a[i][i];
        }


        System.out.println("Gauss");
        for (int i = 0; i < n; i++) {
            System.out.print(x[i] + "\t\t");
        }
        System.out.println("");
    }


    private void Doolittle(double a[][], double b[]) {
        int n = b.length;
        double[][] l = new double[n][n];
        double[][] u = new double[n][n];
        for (int i = 0; i < n; i++) {
            l[i][i] = 1;
        }
        System.arraycopy(a[0], 0, u[0], 0, n);
        for (int i = 1; i < n; i++) {
            l[i][0] = a[i][0] / a[0][0];
        }

        for (int k = 1; k < n; k++) {
            for (int j = k; j < n; j++) {
                double count = 0;
                for (int t = 0; t <= k; t++) {
                    count += l[k][t] * u[t][j];
                }
                u[k][j] = a[k][j] - count;
            }
            for (int i = k + 1; i < n; i++) {
                double count2 = 0;
                for (int t = 0; t <= k - 1; t++) {
                    count2 += l[i][t] * u[t][k];
                }
                l[i][k] = (a[i][k] - count2) / u[k][k];
            }
        }

        double[] y = new double[n];
        double[] x = new double[n];
        y[0] = b[0];
        for (int i = 1; i < n; i++) {
            double count = 0;
            for (int j = 0; j <= i - 1; j++) {
                count += l[i][j] * y[j];
            }
            y[i] = b[i] - count;
        }
        x[n - 1] = y[n - 1] / u[n - 1][n - 1];

        for (int i = n - 2; i >= 0; i--) {
            double count = 0;
            for (int j = i + 1; j <= n - 1; j++) {
                count += u[i][j] * x[j];
            }
            x[i] = (y[i] - count) / u[i][i];
        }
        System.out.println("Doolittle");
        for (int i = 0; i < n; i++) {
            System.out.print(x[i] + "\t\t");
        }
        System.out.println(" ");
    }


    private void jacobi(double a[][], double b[]) {
        int n = b.length;
        double[] Xk = new double[n];
        double[] Xk1 = new double[n];
        for (int i = 0; i < n; i++) {
            Xk[i] = 1;
        }
        int time = 0;
        while (true) {
            time++;
            for (int i = 0; i < n; i++) {
                double count = 0;
                for (int j = 0; j <= i - 1; j++) {
                    count += a[i][j] * Xk[j];
                }
                double count2 = 0;
                for (int j = i + 1; j < n; j++) {
                    count2 += a[i][j] * Xk[j];
                }
                Xk1[i] = (b[i] - count - count2) / a[i][i];
            }
            boolean isStop = true;
            for (int i = 0; i < n; i++) {
                if (Math.abs(Xk1[i] - Xk[i]) > Math.pow(10, -8)) {
                    isStop = false;
                }
            }
            System.arraycopy(Xk1, 0, Xk, 0, n);
            if (isStop) {
                System.out.println("jacobi" + "\t" + "迭代次数：" + time);
                for (int i = 0; i < n; i++) {
                    System.out.print(Xk1[i] + "\t\t");
                }
                System.out.println("");
                break;
            }
        }
    }

    private void GaussSeidel(double a[][], double b[]) {
        int n = b.length;
        double[] Xk = new double[n];
        double[] Xk1 = new double[n];
        for (int i = 0; i < n; i++) {
            Xk[i] = 1;
        }
        int time = 0;
        while (true) {
            time++;
            for (int i = 0; i < n; i++) {
                double count = 0;
                for (int j = 0; j <= i - 1; j++) {
                    count += a[i][j] * Xk1[j];
                }
                double count2 = 0;
                for (int j = i + 1; j < n; j++) {
                    count2 += a[i][j] * Xk[j];
                }
                Xk1[i] = (b[i] - count - count2) / a[i][i];
            }
            boolean isStop = true;
            for (int i = 0; i < n; i++) {
                if (Math.abs(Xk1[i] - Xk[i]) > Math.pow(10, -8)) {
                    isStop = false;
                }
            }
            System.arraycopy(Xk1, 0, Xk, 0, n);
            if (isStop) {
                System.out.println("GaussSeidel" + "\t" + "迭代次数：" + time);
                for (int i = 0; i < n; i++) {
                    System.out.print(Xk1[i] + "\t\t");
                }
                System.out.println("");
                break;
            }
        }
    }


    private void SuccessiveOverRelaxation(double a[][], double b[], double w) {
        int n = b.length;
        double[] Xk = new double[n];
        double[] Xk1 = new double[n];
        for (int i = 0; i < n; i++) {
            Xk[i] = 1;
        }
        int time = 0;
        while (true) {
            time++;
            for (int i = 0; i < n; i++) {
                double count = 0;
                for (int j = 0; j <= i - 1; j++) {
                    count += a[i][j] * Xk1[j];
                }
                double count2 = 0;
                for (int j = i + 1; j < n; j++) {
                    count2 += a[i][j] * Xk[j];
                }
                Xk1[i] = (1 - w) * Xk[i] + w * ((b[i] - count - count2)) / a[i][i];
            }
            boolean isStop = true;
            for (int i = 0; i < n; i++) {
                if (Math.abs(Xk1[i] - Xk[i]) > Math.pow(10, -8)) {
                    isStop = false;
                }
            }
            System.arraycopy(Xk1, 0, Xk, 0, n);
            if (isStop) {
                System.out.println("SuccessiveOverRelaxation" + "\t" + "迭代次数：" + time);
                for (int i = 0; i < n; i++) {
                    System.out.print(Xk1[i] + "\t\t");
                }
                System.out.println("");
                break;
            }
        }
    }

    private  void ConjugateGradient(double A[][], double b[]) {
        int n = b.length;
        double[][] Xk = new double[n][1];
        for (int i = 0; i < n; i++) {
            Xk[i][0] = 0;
        }
        double[][] Xk1 = new double[n][1];
        double[][] B = new double[n][1];
        for (int i = 0; i < n; i++) {
            B[i][0] = b[i];
        }
        Matrix matrixA = new Matrix(A);
        Matrix matrixXk = new Matrix(Xk);
        Matrix matrixXk1 = new Matrix(Xk1);
        Matrix matrixb = new Matrix(B);

        Matrix matrixr = matrixb.minus(matrixA.times(matrixXk));
        Matrix matrixp = matrixr;

        int time = 0;
        while (true) {
            time++;

            double alpha = (matrixr.transpose().times(matrixr).get(0, 0)) / (matrixp.transpose().times(matrixA).times(matrixp).get(0, 0));

            matrixXk1 = matrixXk.plus(matrixp.times(alpha));

            Matrix matrixTempr = matrixr;
            matrixr = matrixr.minus(matrixA.times(matrixp).times(alpha));

            double beta = matrixr.transpose().times(matrixr).get(0, 0) / matrixTempr.transpose().times(matrixTempr).get(0, 0);

            matrixp = matrixr.minus(matrixp.times(beta));

            boolean isStop = true;
            for (int i = 0; i < n; i++) {
                if (Math.abs(matrixXk.get(i, 0) - matrixXk1.get(i, 0)) > Math.pow(10, -8)) {
                    isStop = false;
                }
            }
            if (isStop) {
                System.out.println("ConjugateGradient" + "\t" + "迭代次数：" + time);
                for (int i = 0; i < n; i++) {
                    System.out.print(matrixXk1.get(i, 0) + "\t\t");
                }
                System.out.println("");
                break;
            }
            matrixXk = matrixXk1;
        }
    }
}
