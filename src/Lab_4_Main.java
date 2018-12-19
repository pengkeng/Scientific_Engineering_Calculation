import java.util.function.Function;

public class Lab_4_Main {
    public void start() {

        /**
         * 定义多项式 f = 1 /(1 + x);
         */
        DerivedFunction function = new DerivedFunction(new Function<Double, Double>() {
            @Override
            public Double apply(Double x) {
                return 1f / (1f + x);
            }
        });

        System.out.print("Simpson公式：" + simpson(function, 0, 1, 5));


    }


    private double simpson(DerivedFunction function, double a, double b) {
        if (b < a) {
            double temp = a;
            a = b;
            b = temp;
        }
        return (b - a) / 6f * (function.apply(a) + 4f * function.apply((b - a) / 2f) + function.apply(b));
    }

    private double simpson(DerivedFunction function, double a, double b, int n) {
        if (b < a) {
            double temp = a;
            a = b;
            b = temp;
        }
        double h = (b - a) / n;
        double s = (function.apply(a) - function.apply(b)) / 2;

        for (int j = 1; j <= n; j++) {
            s += 2 * function.apply(a + (j - 0.5) * h) + function.apply(a + j * h);
        }
        return h * s / 3;
    }


    private class DerivedFunction implements Function<Double, Double> {
        private Function<Double, Double> function;

        private DerivedFunction(Function<Double, Double> function) {
            this.function = function;
        }

        @Override
        public Double apply(Double x) {
            return function.apply(x);
        }
    }

}
