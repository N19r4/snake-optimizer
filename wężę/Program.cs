using System;

public interface IOptimizationAlgorithm
{
    string SnakeOptimizer { get; set; }
    double rand()
    {
        Random rand = new Random();
        return rand.NextDouble();
    }
    double Solve()
    {
        // parameters
        int N = 30; // population size
        int dim = 30;
        int T = 500; // max number of iterations
        double lb = -10;
        double ub = 10;

        double fval = 0;

        Func<double[], double> fobj = (x) => Math.Pow(x[0], 2) + Math.Pow(x[1], 2); // replace with objective functions

        // initialization
        int[] vec_flag = { 1, -1 };
        double explorationThreshold = 0.25;
        double exploitationThreshold = 0.6;
        double C1 = 0.5; // constant (given in article)
        double C2 = 0.05; // constant (given in article)
        double C3 = 2; // constant (given in article)

        // create population
        double[][] X = new double[N][];
        Random random = new Random();
        for (int i = 0; i < N; i++)
        {
            X[i] = new double[dim];
            for (int j = 0; j < dim; j++)
            {
                X[i][j] = lb + random.NextDouble() * (ub - lb);
            }
        }

        // calculate objective function values
        double[] fitness = new double[N];
        for (int i = 0; i < N; i++)
        {
            fitness[i] = fobj(X[i]);
        }

        double GYbest = fitness.Min();
        int gbest = Array.IndexOf(fitness, GYbest);
        double[] Xfood = X[gbest];

        // Divide the swarm into two equal groups males and females
        int Nm = N / 2;
        int Nf = N - Nm;
        double[,] Xm = new double[Nm, dim];
        double[,] Xf = new double[Nf, dim];
        double[] fitness_m = new double[Nm];
        double[] fitness_f = new double[Nf];
        Array.Copy(X, 0, Xm, 0, Nm);
        Array.Copy(X, Nm, Xf, 0, Nf);
        Array.Copy(fitness, 0, fitness_m, 0, Nm);
        Array.Copy(fitness, Nm, fitness_f, 0, Nf);

        // find best male
        double fitnessBest_m = fitness_m.Min();
        int gbest1 = Array.IndexOf(fitness_m, fitnessBest_m);
        double[] Xbest_m = Enumerable.Range(0, Xm.GetLength(1))
                .Select(x => Xm[gbest1, x])
                .ToArray();


        // find best feamle
        double fitnessBest_f = fitness_f.Min();
        int gbest2 = Array.IndexOf(fitness_f, fitnessBest_f);
        double[] Xbest_f = Enumerable.Range(0, Xf.GetLength(1))
                .Select(x => Xf[gbest2, x])
                .ToArray();

        // start iterating
        for (int t = 1; t <= T; t++)
        {
            double Temp = Math.Exp(-((double)t / T));  //eq.(4)
            double Q = C1 * Math.Exp(((double)(t - T) / T)); //eq.(5)
            if (Q > 1)
                Q = 1;

            double[,] Xnewm = new double[Nm, dim];
            double[,] Xnewf = new double[Nf, dim];

            

            // Exploration Phase (no Food)
            if (Q < explorationThreshold)
            {
                for (int i = 0; i < Nm; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        int rand_leader_index = (int)(Nm * new Random().NextDouble() + 1);
                        double[] X_randm = Enumerable.Range(0, Xm.GetLength(1))
                .Select(x => Xm[rand_leader_index, x])
                .ToArray(); ;
                        int flag_index = (int)(2 * new Random().NextDouble() + 1);
                        double Flag = vec_flag[flag_index];
                        double Am = Math.Exp(-fitness_m[rand_leader_index] / (fitness_m[i] + double.Epsilon)); //eq.(7)
                        Xnewm[i, j] = X_randm[j] + Flag * C2 * Am * ((ub - lb) * new Random().NextDouble() + lb); //eq.(6)
                    }
                }

                for (int i = 0; i < Nf; i++)
                {
                    for (int j = 0; j < dim; j++)
                    {
                        int rand_leader_index = (int)(Nf * new Random().NextDouble() + 1);
                        double[] X_randf = Enumerable.Range(0, Xf.GetLength(1))
                .Select(x => Xf[rand_leader_index, x])
                .ToArray(); ;
                        int flag_index = (int)(2 * new Random().NextDouble() + 1);
                        double Flag = vec_flag[flag_index];
                        double Af = Math.Exp(-fitness_f[rand_leader_index] / (fitness_f[i] + double.Epsilon)); //eq.(9)
                        Xnewf[i, j] = X_randf[j] + Flag * C2 * Af * ((ub - lb) * new Random().NextDouble() + lb); //eq.(8)
                    }
                }
            }
            else //Exploitation Phase (Food Exists)
            {
                if (Temp > exploitationThreshold) //hot
                {
                    for (int i = 0; i < Nm; i++)
                    {
                        int flag_index = (int)(2 * rand() + 1);
                        int Flag = vec_flag[flag_index];
                        for (int j = 0; j < dim; j++)
                        {
                            Xnewm[i, j] = Xfood[j] + C3 * Flag * Temp * rand() * (Xfood[j] - Xm[i, j]); //eq.(10)
                        }
                    }
                    for (int i = 0; i < Nf; i++)
                    {
                        int flag_index = (int)(2 * rand() + 1);
                        int Flag = vec_flag[flag_index];
                        for (int j = 0; j < dim; j++)
                        {
                            Xnewf[i, j] = Xfood[j] + Flag * C3 * Temp * rand() * (Xfood[j] - Xf[i, j]); //eq.(10)
                        }
                    }
                }
                else //cold
                {
                    if (rand() > 0.6) //fight mode
                    {
                        for (int i = 0; i < Nm; i++)
                        {
                            for (int j = 0; j < dim; j++)
                            {
                                double FM = Math.Exp(-(fitnessBest_f) / (fitness_m[i] + double.Epsilon)); //eq.(13)
                                Xnewm[i, j] = Xm[i, j] + C3 * FM * rand() * (Q * Xbest_f[j] - Xm[i, j]); //eq.(11)
                            }
                        }
                        for (int i = 0; i < Nf; i++)
                        {
                            for (int j = 0; j < dim; j++)
                            {
                                double FF = Math.Exp(-(fitnessBest_m) / (fitness_f[i] + double.Epsilon)); //eq.(14)
                                Xnewf[i, j] = Xf[i, j] + C3 * FF * rand() * (Q * Xbest_m[j] - Xf[i, j]); //eq.(12)
                            }
                        }
                    }
                    else //mating mode
                    {
                        for (int i = 0; i < Nm; i++)
                        {
                            for (int j = 0; j < dim; j++)
                            {
                                double Mm = Math.Exp(-fitness_f[i] / (fitness_m[i] + double.Epsilon)); //eq.(17)
                                Xnewm[i, j] = Xm[i, j] + C3 * rand() * Mm * (Q * Xf[i, j] - Xm[i, j]); //eq.(15)
                            }
                        }
                        for (int i = 0; i < Nf; i++)
                        {
                            for (int j = 0; j < dim; j++)
                            {
                                double Mf = Math.Exp(-fitness_m[i] / (fitness_f[i] + double.Epsilon)); // eq.(18)
                                Xnewf[i, j] = Xf[i, j] + C3 * rand() * Mf * (Q * Xm[i, j] - Xf[i, j]); // eq.(16)
                            }
                        }

                        int flag_index = (int)Math.Floor(2 * rand() + 1);
                        int egg = vec_flag[flag_index];
                        
                        // if egg hatch, replace worst 
                        if (egg == 1)
                        {
                            int gworst;
                            double GYworst = fitness_m.Max();
                            gworst = Array.IndexOf(fitness_m, GYworst);
                            for (int j = 0; j < dim; j++)
                            {
                                Xnewm[gworst, j] = lb + rand() * (ub - lb); // eq.(19)
                            }
                            GYworst = fitness_f.Max();
                            gworst = Array.IndexOf(fitness_f, GYworst);
                            for (int j = 0; j < dim; j++)
                            {
                                Xnewf[gworst, j] = lb + rand() * (ub - lb); // eq.(20)
                            }
                        }
                    }
                }
            }

            for (int j = 0; j < Nm; j++)
            {
                bool[] Flag4ub = new bool[Xnewm.GetLength(1)];
                bool[] Flag4lb = new bool[Xnewm.GetLength(1)];
                for (int i = 0; i < Xnewm.GetLength(1); i++)
                {
                    Flag4ub[i] = Xnewm[j, i] > ub[i];
                    Flag4lb[i] = Xnewm[j, i] < lb[i];
                }
                for (int i = 0; i < Xnewm.GetLength(1); i++)
                {
                    Xnewm[j, i] = (Xnewm[j, i] * (!(Flag4ub[i] || Flag4lb[i]))) + ub[i] * Convert.ToInt32(Flag4ub[i]) + lb[i] * Convert.ToInt32(Flag4lb[i]);
                }
                double y = feval(fobj, Xnewm.GetRow(j));
                if (y < fitness_m[j])
                {
                    fitness_m[j] = y;
                    Xm[j] = Xnewm.GetRow(j);
                }
            }

            [Ybest2, gbest2] = min(fitness_f);

            if (Ybest1 < fitnessBest_m)
            {
                Xbest_m = Xm(gbest1,:);
                fitnessBest_m = Ybest1;
            }
            if (Ybest2 < fitnessBest_f)
            {
                Xbest_f = Xf(gbest2,:);
                fitnessBest_f = Ybest2;
            }

            if (Ybest1 < Ybest2)
            {
                gbest_t(t) = min(Ybest1);
            }
            else
            {
                gbest_t(t) = min(Ybest2);
            }

            if (fitnessBest_m < fitnessBest_f)
            {
                GYbest = fitnessBest_m;
                Xfood = Xbest_m;
            }
            else
            {
                GYbest = fitnessBest_f;
                Xfood = Xbest_f;
            }
            fval = GYbest;

        }

        return fval;
    }
    double[] XBest { get; set; }
    double FBest { get; set; }
    int NumberOfEvaluationFitnessFunction { get; set; }
 }



class TestClass : IOptimizationAlgorithm
{
    public string SnakeOptimizer { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
    public double[] XBest { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
    public double FBest { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }
    public int NumberOfEvaluationFitnessFunction { get => throw new NotImplementedException(); set => throw new NotImplementedException(); }

    static void Main(string[] args)
    {
        
    }
}
