#include <vector>
#include <cmath>
extern "C" std::vector<double> jcalc(std::vector<double> listvalue, int curstate)
{
    std::vector<double> ret;
    if (curstate == 1)
    {
        double xp= listvalue[0];
        double yp= listvalue[1];
        double xd= listvalue[2];
        double yd= listvalue[3];
        double Entry = 0;
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=1;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=1;
        ret.push_back(Entry);
        Entry=-0.0575997658817729;
        ret.push_back(Entry);
        Entry=0.000200959896519766;
        ret.push_back(Entry);
        Entry=-2.89995083970656;
        ret.push_back(Entry);
        Entry=0.00877200894463775;
        ret.push_back(Entry);
        Entry=-0.000174031357370456;
        ret.push_back(Entry);
        Entry=-0.0665123984901026;
        ret.push_back(Entry);
        Entry=-0.00875351105536225;
        ret.push_back(Entry);
        Entry=-2.90300269286856;
        ret.push_back(Entry);
        return ret;
    }
    if (curstate == 2)
    {
        double xp= listvalue[0];
        double yp= listvalue[1];
        double xd= listvalue[2];
        double yd= listvalue[3];
        double Entry = 0;
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=1;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=1;
        ret.push_back(Entry);
        Entry=-0.575999943070835;
        ret.push_back(Entry);
        Entry=0.000262486079431672;
        ret.push_back(Entry);
        Entry=-19.2299795908647;
        ret.push_back(Entry);
        Entry=0.00876275931760007;
        ret.push_back(Entry);
        Entry=-0.000262486080737868;
        ret.push_back(Entry);
        Entry=-0.575999940191886;
        ret.push_back(Entry);
        Entry=-0.00876276068239993;
        ret.push_back(Entry);
        Entry=-19.2299765959399;
        ret.push_back(Entry);
        return ret;
    }
    if (curstate == 3)
    {
        double xp= listvalue[0];
        double yp= listvalue[1];
        double xd= listvalue[2];
        double yd= listvalue[3];
        double Entry = 0;
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=1;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=1;
        ret.push_back(Entry);
        Entry=5.75894721132000e-5;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0.00876276000000000;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        Entry=-0.00876276000000000;
        ret.push_back(Entry);
        Entry=0;
        ret.push_back(Entry);
        return ret;
    }

}

extern "C" void thinHandle(int state, std::vector<double> &notbloating, std::vector<double> &thin_indicate, std::vector<double> &bloating, std::vector<double> &not_thin)
{
    if (int(state) == 1)
    {
        notbloating=std::vector<double>{4,5};
        thin_indicate=std::vector<double>{};
        bloating=std::vector<double>{0,1,2,3};
        not_thin=std::vector<double>{0,1,2,3};
        return;
    }
    if (int(state) == 2)
    {
        notbloating=std::vector<double>{4,5};
        thin_indicate=std::vector<double>{};
        bloating=std::vector<double>{0,1,2,3};
        not_thin=std::vector<double>{0,1,2,3};
        return;
    }
    if (int(state) == 3)
    {
        notbloating=std::vector<double>{4,5};
        thin_indicate=std::vector<double>{};
        bloating=std::vector<double>{0,1,2,3};
        not_thin=std::vector<double>{0,1,2,3};
        return;
    }
}
