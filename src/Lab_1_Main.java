

public class Lab_1_Main {

    /**
     * 二分法 牛顿法 牛顿下山法 弦截法
     * x^2 - 3 * x + 2 - e^x = 0;
     * x^3 - x - 1 = 0;
     * 当区间逼近与 10^-8 计算结束
     */

    public void start() {
        dichotomy1();
        newton1();
        newtonDown1();
        secant1();
        dichotomy2();
        newton2();
        newtonDown2();
        secant2();
        double[][] a= new double[1][1];
    }

    /**
     * x^2 - 3 * x + 2 - e^x = 0
     *
     * @param x
     * @return
     */
    private double fun1(double x) {
        return Math.pow(x, 2) - 3 * x + 2 - Math.pow(Math.E, x);
    }

    /**
     * 2 * x - 3 - e^x;
     *
     * @param x
     * @return
     */
    private double derivative_fun1(double x) {
        return 2 * x - 3 - Math.pow(Math.E, x);
    }

    /**
     * x^3 - x - 1 = 0
     *
     * @param x
     * @return
     */
    private double fun2(double x) {
        return Math.pow(x, 3) - x - 1;
    }

    /**
     * 3 * x - 1
     *
     * @param x
     * @return
     */
    private double derivative_fun2(double x) {
        return 3 * Math.pow(x, 2) - 1;
    }

    /**
     * 二分法
     */
    private void dichotomy1() {
        //已知 fun(0) = 1,fun(1) = -e;
        //因此初始 x1 = 0; x2 =1;
        double x1 = 0, x2 = 1;
        double x3 = (x1 + x2) / 2;
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x3:" + x3 + " 函数值 ：" + fun1(x3));
            x3 = (x1 + x2) / 2;
            if (Math.abs(x1 - x2) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
            if (fun1(x3) < 0) {
                x2 = x3;
            } else if (fun1(x3) > 0) {
                x1 = x3;
            } else {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
            }
        }
    }

    /**
     * 二分法
     */
    private void dichotomy2() {
        //已知 fun1(1) < 0,fun1(2) > 0;
        //因此初始 x1 = 0; x2 =1;
        double x1 = 1, x2 = 2;
        double x3 = (x1 + x2) / 2;
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x3:" + x3 + " 函数值 ：" + fun2(x3));
            x3 = (x1 + x2) / 2;
            if (Math.abs(x1 - x2) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
            if (fun2(x3) > 0) {
                x2 = x3;
            } else if (fun2(x3) < 0) {
                x1 = x3;
            } else {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
            }
        }
    }

    /**
     * 牛顿收敛法
     */
    private void newton1() {
        // 去x0初值为 0；
        // x0 = 0;
        double x1 = 0;
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x1:" + x1 + " 函数值 ：" + fun1(x1));
            double x0 = x1;
            x1 = x1 - fun1(x1) / derivative_fun1(x1);
            if (Math.abs(x1 - x0) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x1 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
        }
    }

    /**
     * 牛顿收敛法
     */
    private void newton2() {
        // 去x0初值为 0；
        // x0 = 0;
        double x1 = 0;
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x1:" + x1 + " 函数值 ：" + fun2(x1));
            double x0 = x1;
            x1 = x1 - fun2(x1) / derivative_fun2(x1);
            if (Math.abs(x1 - x0) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x1 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
        }
    }

    /**
     * 牛顿下山法
     */
    private void newtonDown1() {
        // 去x0初值为 0；
        // x0 = 0;
        double x1 = 0;
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x1:" + x1 + " 函数值 ：" + fun1(x1));
            double x0 = x1;
            x1 = x1 - fun1(x1) / derivative_fun1(x1);
            if (Math.abs(fun1(x0)) < (Math.abs(fun1(x1)))) {
                x1 = x0 - 0.5 * fun1(x0) / derivative_fun1(x0);
            }
            if (Math.abs(x1 - x0) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x1 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
        }
    }

    /**
     * 牛顿下山法
     */
    private void newtonDown2() {
        // 去x0初值为 0；
        // x0 = 0;
        double x1 = 0;
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x1:" + x1 + " 函数值 ：" + fun2(x1));
            double x0 = x1;
            x1 = x1 - fun2(x1) / derivative_fun2(x1);
            if (Math.abs(fun2(x0)) < (Math.abs(fun2(x1)))) {
                x1 = x0 - 0.5 * fun2(x0) / derivative_fun2(x0);
            }
            if (Math.abs(x1 - x0) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x1 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
        }
    }

    /**
     * fun(0) = 1,fun(1) = -e;
     * x1 = 0; x2 =1;
     */
    private void secant1() {
        //已知 fun(0) = 1,fun(1) = -e;
        //因此初始 x1 = 0; x2 =1;
        double x1 = 0, x2 = 1;
        int time = 0;
        double x3 = x1 - fun1(x1) / (fun1(x2) - fun1(x1)) * (x2 - x1);
        ;
        while (true) {
            time ++;
            System.out.println("x3:" + x3 + " 函数值 ：" + fun1(x3));
            x3 = x1 - fun1(x1) / (fun1(x2) - fun1(x1)) * (x2 - x1);
            if (Math.abs(x1 - x2) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
            if (fun1(x3) < 0) {
                x2 = x3;
            } else if (fun1(x3) > 0) {
                x1 = x3;
            } else {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
        }
    }

    /**
     * fun1(1) < 0,fun1(2) > 0;
     * x1 = 0; x2 =1;
     */
    private void secant2() {
        double x1 = 1, x2 = 2;
        double x3 = x1 - fun2(x1) / (fun2(x2) - fun2(x1)) * (x2 - x1);
        int time = 0;
        while (true) {
            time ++;
            System.out.println("x3:" + x3 + " 函数值 ：" + fun2(x3));
            x3 = x1 - fun2(x1) / (fun2(x2) - fun2(x1)) * (x2 - x1);
            if (Math.abs(x1 - x2) < (Math.pow(10, -8))) {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
            if (fun2(x3) > 0) {
                x2 = x3;
            } else if (fun2(x3) < 0) {
                x1 = x3;
            } else {
                System.out.println("<---------------------------------------------->");
                System.out.println("方程式的实根近似为：" + x3 + " 迭代了 ：" + time +" 次");
                System.out.println("<---------------------------------------------->");
                break;
            }
        }
    }
}
