namespace Pricer.Numerics;

public static class BinomialTree
{
    public static double EuropeanOptionPrice(
        OptionType optionType,
        double underlyingPrice,
        double strike,
        double riskFreeRate,
        double timeToMaturity,
        double volatility,
        int steps)
    {
        if (underlyingPrice <= 0.0)
            throw new ArgumentException("Underlying price must be positive.");

        if (strike <= 0.0)
            throw new ArgumentException("Strike must be positive.");

        if (timeToMaturity < 0.0)
            throw new ArgumentException("Time to maturity cannot be negative.");

        if (volatility < 0.0)
            throw new ArgumentException("Volatility cannot be negative.");

        if (steps <= 0)
            throw new ArgumentException("Number of steps must be positive.");

        if (timeToMaturity == 0.0)
        {
            return optionType == OptionType.Call
                ? Math.Max(underlyingPrice - strike, 0.0)
                : Math.Max(strike - underlyingPrice, 0.0);
        }

        if (volatility == 0.0)
        {
            double deterministicTerminalStockPrice = underlyingPrice * Math.Exp(riskFreeRate * timeToMaturity);

            double terminalPayoff = optionType == OptionType.Call
                ? Math.Max(deterministicTerminalStockPrice - strike, 0.0)
                : Math.Max(strike - deterministicTerminalStockPrice, 0.0);

            return Math.Exp(-riskFreeRate * timeToMaturity) * terminalPayoff;
        }

        double dt = timeToMaturity / steps;
        double discount = Math.Exp(-riskFreeRate * dt);

        double u = Math.Exp(volatility * Math.Sqrt(dt));
        double d = Math.Exp(-volatility * Math.Sqrt(dt));

        double q = (Math.Exp(riskFreeRate * dt) - d) / (u - d);

        if (q < 0.0 || q > 1.0)
            throw new InvalidOperationException("Risk-neutral probability is outside [0, 1]. Check the input parameters.");

        double[] optionValues = new double[steps + 1];

        // Terminal payoffs at maturity T.
        // Node j means: j up moves and steps - j down moves.
        for (int j = 0; j <= steps; j++)
        {
            double stockPrice = underlyingPrice * Math.Pow(u, j) * Math.Pow(d, steps - j);

            optionValues[j] = optionType == OptionType.Call
                ? Math.Max(stockPrice - strike, 0.0)
                : Math.Max(strike - stockPrice, 0.0);
        }

        // Backward induction.
        for (int i = steps - 1; i >= 0; i--)
        {
            for (int j = 0; j <= i; j++)
            {
                optionValues[j] = discount * (q * optionValues[j + 1] + (1.0 - q) * optionValues[j]);
            }
        }

        return optionValues[0];
    }

    public static double CRRUpFactor(double volatility, double timeToMaturity, int steps)
    {
        ValidateTreeParameters(volatility, timeToMaturity, steps);
        double dt = timeToMaturity / steps;
        return Math.Exp(volatility * Math.Sqrt(dt));
    }

    public static double CRRDownFactor(double volatility, double timeToMaturity, int steps)
    {
        ValidateTreeParameters(volatility, timeToMaturity, steps);
        double dt = timeToMaturity / steps;
        return Math.Exp(-volatility * Math.Sqrt(dt));
    }

    public static double RiskNeutralProbability(double riskFreeRate, double volatility, double timeToMaturity, int steps)
    {
        ValidateTreeParameters(volatility, timeToMaturity, steps);

        if (volatility == 0.0)
            throw new InvalidOperationException("Risk-neutral probability is not well-defined when volatility is zero.");

        double dt = timeToMaturity / steps;
        double u = Math.Exp(volatility * Math.Sqrt(dt));
        double d = Math.Exp(-volatility * Math.Sqrt(dt));

        return (Math.Exp(riskFreeRate * dt) - d) / (u - d);
    }

    private static void ValidateTreeParameters(double volatility, double timeToMaturity, int steps)
    {
        if (volatility < 0.0)
            throw new ArgumentException("Volatility cannot be negative.");

        if (timeToMaturity <= 0.0)
            throw new ArgumentException("Time to maturity must be positive.");

        if (steps <= 0)
            throw new ArgumentException("Number of steps must be positive.");
    }
}
