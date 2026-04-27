using MathNet.Numerics.Distributions; // gebruik een externe bibliotheek om met kansverdelingen te werken


namespace Pricer.Numerics; // mapnaam

public enum OptionType // lijst van vaste keuzes
{
    Call,
    Put
}

public class BlackScholes
{   // parameters voor het BS model, private = alleen binnen deze class te gebruiken, double = kommagetal
    private OptionType option;
    private double r; // Risk-free interest rate
    private double T; // Time to maturity
    private double sigma; // Volatility
    private double K; // Strike price
    private double S; // Underlying asset price
    private double q; // Dividend yield

    // constructor = functie die wordt uitgevoerd wanneer we een nieuw BS-object maken
    public BlackScholes(
        OptionType optionType,
        double riskFreeRate,
        double timeToMaturity,
        double volatility,
        double strike,
        double underlyingPrice,
        double dividendYield = 0.0) // standaardwaarde, als je geen dividend geeft wordt automatisch 0 gebruikt
    {
        option = optionType;
        r = riskFreeRate;
        T = timeToMaturity;
        sigma = volatility;
        K = strike;
        S = underlyingPrice;
        q = dividendYield;
    }

    // Cumulative normal distribution
    private static double N(double x) // static = algemene hulpfunctie, gebruikt de pm niet
        => Normal.CDF(0.0, 1.0, x); // gemiddelde = 0, standaarddeviatie = 1, x = waarde waarvoor we de cumulatieve kans willen weten

    // Normal probability density function
    private static double n(double x)
        => Normal.PDF(0.0, 1.0, x);


    // functie op de Black–Scholes prijs van een optie te berekenen
    public double Price() // public = deze functie kan van buiten de class worden aangeroepen
    {
        // randgevallen
        if (T <= 0.0) // maturity nul of negatief = optie vervalt nu meteen
        {
            return option == OptionType.Call // check of option = call
                ? Math.Max(S - K, 0.0) // ja geef dit
                : Math.Max(K - S, 0.0); // nee het is een put option dus geef dit
        }

        if (sigma <= 0.0) // volatiliteit nul of negatief = geen onzekerheid
        {
            double forwardIntrinsic = option == OptionType.Call
                ? Math.Max(S * Math.Exp(-q * T) - K * Math.Exp(-r * T), 0.0)
                : Math.Max(K * Math.Exp(-r * T) - S * Math.Exp(-q * T), 0.0);

            return forwardIntrinsic;
        }

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
        double d2 = d1 - sigma * sqrtT;

        if (option == OptionType.Call) // als het een Call optie is, doe dan dit
        {
            return S * Math.Exp(-q * T) * N(d1)
                 - K * Math.Exp(-r * T) * N(d2);
        }
        else // als het een Put optie is, doe dan dit
        {
            return K * Math.Exp(-r * T) * N(-d2)
                 - S * Math.Exp(-q * T) * N(-d1);
        }
    }



    // Delta: derivative of option price with respect to the underlying price S
    public double Delta()
    {
        if (T <= 0.0)
        {
            if (option == OptionType.Call)
                return S > K ? 1.0 : 0.0;

            return S < K ? -1.0 : 0.0;
        }

        if (sigma <= 0.0)
            return option == OptionType.Call
                ? (S * Math.Exp((r - q) * T) > K ? Math.Exp(-q * T) : 0.0)
                : (S * Math.Exp((r - q) * T) < K ? -Math.Exp(-q * T) : 0.0);

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);

        return option == OptionType.Call
            ? Math.Exp(-q * T) * N(d1)
            : Math.Exp(-q * T) * (N(d1) - 1.0);
    }

    // Gamma: derivative of delta with respect to the underlying price S
    public double Gamma()
    {
        if (T <= 0.0 || sigma <= 0.0)
            return 0.0;

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);

        return Math.Exp(-q * T) * n(d1) / (S * sigma * sqrtT);
    }

    // Theta: derivative of option price with respect to calendar time.
    // This is the standard market convention: theta is usually negative for long options.
    public double Theta()
    {
        if (T <= 0.0 || sigma <= 0.0)
            return 0.0;

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
        double d2 = d1 - sigma * sqrtT;

        double firstTerm = -S * Math.Exp(-q * T) * n(d1) * sigma / (2.0 * sqrtT);

        if (option == OptionType.Call)
        {
            return firstTerm
                 - r * K * Math.Exp(-r * T) * N(d2)
                 + q * S * Math.Exp(-q * T) * N(d1);
        }

        return firstTerm
             + r * K * Math.Exp(-r * T) * N(-d2)
             - q * S * Math.Exp(-q * T) * N(-d1);
    }

    // Derivative of price with respect to volatility
    public double Vega() // de vega van de optie berekenen (vega = hoe gevoelig is de prijs van de optie voor veranderingen in de volatiliteit)
    {
        if (T <= 0.0 || sigma <= 0.0) // || = of
            return 0.0;

        double sqrtT = Math.Sqrt(T);
        double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);

        return S * Math.Exp(-q * T) * sqrtT * n(d1); // formule voor vega
    }
}
    // Test price compared with party at the mooonlight option calculator
    // Test put call parity for consistency
    // --> zie program.cs


    public static class ImpliedVolatilityCalculator
    {
        public static double Compute(
            OptionType optionType,
            double marketPrice, // P_mkt
            double interestRate,
            double timeToMaturity,
            double strike,
            double underlyingPrice,
            double initialGuess = 0.2, // startgok voor de volatitiliteit
            double tolerance = 1e-8, // wanneer stoppen we (epsilon)
            int maxIterations = 100)
        {
            double sigma = initialGuess; // variabele sigma aanmaken en gelijkstellen aan de startgok (huidige schatting van de vol)

            for (int i = 0; i < maxIterations; i++) // herhaal de Newton-stap maximaal maxIterations keer
                                                    // int i = 0 : teller begint op 0, i < maxIterations: doorgaan zolang we onder de limiet zitten, i++: verhoog de teller telkens met 1
            {   // BS object aanmaken met de huidige sigma = BS(sigma)
                var bs = new BlackScholes(
                    optionType,
                    interestRate,
                    timeToMaturity,
                    sigma,
                    strike,
                    underlyingPrice);

                double modelPrice = bs.Price(); // BS(sigma)
                double vega = bs.Vega(); // P_mkt
                double error = modelPrice - marketPrice;

                if (Math.Abs(error) < tolerance)
                    return sigma; // oplossing is goed genoeg, dus geef sigma terug

                if (Math.Abs(vega) < 1e-12) // bescherming tegen te kleine Vega
                    throw new InvalidOperationException("Vega is too small; Newton-Raphson becomes unstable.");

                sigma = sigma - error / vega; // Newton update

                if (sigma <= 0.0) // bescherming tegen negatieve volatiliteit
                    sigma = 1e-8;
            }
            // als de loop niet zou convergeren:
            throw new InvalidOperationException("Black-Scholes implied volatility did not converge.");
        } 
        // testje zie program.cs

        public static double BachelierImpliedVolATM(
            double optionPrice, // C_0^B
            double underlyingPrice, // S(0)
            double timeToMaturity,
            double interestRate)
        {
            _ = interestRate; // "ja ik weet dat deze pm bestaat maar ik gebruik hem bewust niet"

            if (optionPrice <= 0.0)
                throw new ArgumentException("Option price must be strictly positive.");

            if (underlyingPrice <= 0.0)
                throw new ArgumentException("Underlying price must be strictly positive.");

            if (timeToMaturity <= 0.0)
                throw new ArgumentException("Time to maturity must be strictly positive.");

            double sqrtT = Math.Sqrt(timeToMaturity);

            return optionPrice * Math.Sqrt(2.0 * Math.PI) / (underlyingPrice * sqrtT);
        }
    }
