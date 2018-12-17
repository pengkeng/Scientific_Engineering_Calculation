public class Lab_3_Main {

    private LagrangeInterpolationHelper lagrangeInterpolationHelper;
    private NewtonInterpolationHelper newtonInterpolationHelper;


    /**
     * 初始化拉格朗日插值法
     *
     * @param fun
     * @param point
     */
    public void initLagrangeInterpolation(double[] fun, double[] point) {
        int n = fun.length;
        lagrangeInterpolationHelper = new LagrangeInterpolationHelper() {
            @Override
            public double getValue(double x) {
                double value = 0;
                for (int i = 0; i < n; i++) {
                    double molecule = 1;
                    for (int j = 0; j < n; j++) {
                        if (j != i) {
                            molecule *= (x - point[j]);
                        }
                    }
                    double denominator = 1;
                    for (int j = 0; j < n; j++) {
                        if (j != i) {
                            denominator *= (point[i] - point[j]);
                        }
                    }
                    value += molecule / denominator * fun[i];
                }
                return value;
            }
        };
    }


    /**
     * 初始化牛顿插值法
     *
     * @param fun
     * @param point
     */
    public void initNewtonInterpolation(double[] fun, double[] point) {
        int n = fun.length;
        newtonInterpolationHelper = new NewtonInterpolationHelper() {
            @Override
            public double getValue(double x) {
                double value = fun[0];
                for (int i = 1; i < n; i++) {
                    double temp = 1;
                    double[] tempArray = new double[i + 1];

                    double[] tempArrayValue = new double[i + 1];

                    System.arraycopy(point, 0, tempArray, 0, i + 1);
                    System.arraycopy(fun, 0, tempArrayValue, 0, i + 1);
                    temp = difference_Quotient(tempArray, tempArrayValue);
                    for (int j = 0; j < i; j++) {
                        temp *= (x - point[j]);
                    }
                    value += temp;
                }
                return value;
            }
        };
    }

    public LagrangeInterpolationHelper getLagrangeInterpolationHelper() {
        return lagrangeInterpolationHelper;
    }

    public NewtonInterpolationHelper getNewtonInterpolationHelper() {
        return newtonInterpolationHelper;
    }

    public void start() {
        double[] fun = new double[]{0.98, 0.92, 0.81, 0.64, 0.36};
        double[] point = new double[]{0.2, 0.4, 0.6, 0.8, 1.0};
        initLagrangeInterpolation(fun, point);
        System.out.println(getLagrangeInterpolationHelper().getValue(0.2 + 0.08 * 0));
        System.out.println(getLagrangeInterpolationHelper().getValue(0.2 + 0.08 * 1));
        System.out.println(getLagrangeInterpolationHelper().getValue(0.2 + 0.08 * 11));
        System.out.println(getLagrangeInterpolationHelper().getValue(0.2 + 0.08 * 10));
        initNewtonInterpolation(fun, point);
        System.out.println(getNewtonInterpolationHelper().getValue(0.2 + 0.08 * 0));
        System.out.println(getNewtonInterpolationHelper().getValue(0.2 + 0.08 * 1));
        System.out.println(getNewtonInterpolationHelper().getValue(0.2 + 0.08 * 11));
        System.out.println(getNewtonInterpolationHelper().getValue(0.2 + 0.08 * 10));

        double n = 10;
        double[] x = new double[(int) n];
        double[] values = new double[(int) n];
        for (int i = 0; i < n; i++) {
            x[i] = -1 + 2 * i / n;
            values[i] = func(x[i]);
        }
        initLagrangeInterpolation(values, x);
        initNewtonInterpolation(values, x);
        for (int i = 0; i < 20f; i++) {
            double a = -1 + 2 * i / 20f;
            System.out.println(+getLagrangeInterpolationHelper().getValue(a) + "   \t\t" + func(a) + "\t\t" + getNewtonInterpolationHelper().getValue(a));
        }


    }

    private double func(double x) {
        return -(1 / (1 + 25 * x * x));
    }


    public interface LagrangeInterpolationHelper {
        double getValue(double x);
    }

    public interface NewtonInterpolationHelper {
        double getValue(double x);
    }


    /**
     * 差商计算
     *
     * @param x
     * @param value
     * @return
     */
    private double difference_Quotient(double[] x, double[] value) {
        if (x.length == 2) {
            return (value[1] - value[0]) / (x[1] - x[0]);
        } else {

            double[] startx = new double[x.length - 1];
            System.arraycopy(x, 0, startx, 0, x.length - 1);

            double[] startValue = new double[x.length - 1];
            System.arraycopy(value, 0, startValue, 0, x.length - 1);

            double[] endx = new double[x.length - 1];
            System.arraycopy(x, 1, endx, 0, x.length - 1);

            double[] endValue = new double[x.length - 1];
            System.arraycopy(value, 1, endValue, 0, x.length - 1);

            return (difference_Quotient(endx, endValue) - difference_Quotient(startx, startValue)) / (x[x.length - 1] - x[0]);
        }

    }
}
