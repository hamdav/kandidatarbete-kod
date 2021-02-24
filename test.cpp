#include <set>
#include <vector>
#include <iostream>
#include <cassert>
#include <iterator>

struct GM {
    int sup;
    bool sub;
};

bool operator<(GM w1, GM w2) { return w1.sup == w2.sup ? w1.sub < w2.sub : w1.sup < w2.sup; }
bool operator>(GM w1, GM w2) { return w1.sup == w2.sup ? w1.sub > w2.sub : w1.sup > w2.sup; }
bool operator==(GM w1, GM w2) { return w1.sup == w2.sup && w1.sub == w2.sub; }

struct Term {
    int factor;
    std::set<GM> ws;
};

bool operator<(const Term t1, const Term t2) { return t1.ws < t2.ws; }
bool operator>(const Term t1, const Term t2) { return t1.ws > t2.ws; }
Term operator*(const Term t1, const Term t2) { 

    // Merge the two sorted lists t1.ws and t2.ws
    // keeping track of the sign
    unsigned long remainingElemsIn1{ t1.ws.size() };
    bool parity{ 0 };
    std::set<GM> ws{};
    auto it1{ t1.ws.begin() };
    auto it2{ t2.ws.begin() };

    while (it1 != t1.ws.end() && it2 != t2.ws.end()) {
        if (*it1 == *it2) 
            return {0, {}};
        if (*it1 < *it2){
            ws.insert(ws.end(), *it1);
            --remainingElemsIn1;
            ++it1;
        }
        else {
            ws.insert(ws.end(), *it2);
            ++it2;
            parity = (parity + remainingElemsIn1) % 2;
        }
    }
    while (it1 != t1.ws.end()){
        ws.insert(ws.end(), *it1);
        ++it1;
    }
    while (it2 != t2.ws.end()){
        ws.insert(ws.end(), *it2);
        ++it2;
    }

    return {(parity ? -1 : 1) * t1.factor * t2.factor, ws};
}
Term operator+(const Term t1, const Term t2) { 
    assert(t1.ws == t2.ws && "Can't add different terms");
    return {t1.factor + t2.factor, t1.ws};
}
Term operator*(Term t, const GM w) {

    // Try to insert w
    auto itPair{ t.ws.insert(w) };

    // If insertion was successful
    if (itPair.second) {
        bool parity = (std::distance(itPair.first, t.ws.end()) + 1) % 2;
        if (parity)
            t.factor *= -1;
        return t;
    }
    // Otherwise, w was already in t
    else {
        return {0, {}};
    }
}

class Poly{
public:
    std::set<Term> terms;
    Poly(){
        terms = std::set<Term>{Term{1, {}}};
    }
};

std::ostream& operator<<(std::ostream& out, GM w) {
    if (w.sub)
        out << "V";
    else
        out << "U";
    out << "[" << w.sup << "]";
    return out;
}

std::ostream& operator<<(std::ostream& out, Term t) {
    if (t.factor >= 0)
        out << " +" << t.factor;
    else
        out << " -" << -t.factor;
    for (auto w : t.ws) {
        out << w;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, Poly p) {
    out << "(";
    for (auto t : p.terms) 
        out << t;
    out << ")";
    return out;
}


// S is U^a V^a
// Multiply by S is thus for each term adding a term with an aditional factor U^a V^a
// (if there is not already a U^a or V^a)

void multByS(Poly& p, int n) {
    // Create a set for the new terms
    std::set<Term> newTerms{};

    // For each old term, add the relevant new resulting terms. 
    for (auto t : p.terms) {
        for (int i{ 1 }; i <= n; i++){
            Term newt{ t * GM{i, 0} * GM{i, 1} };
            auto foundIt{ newTerms.find(newt) };
            // If it was already there, add the factors
            if (foundIt != newTerms.end()) {
                // If the factors cancel, remove it from the terms
                if (foundIt->factor + t.factor == 0) 
                    newTerms.erase(foundIt);
                // Otherwise add the factors. 
                else {
                    newt.factor += foundIt->factor;
                    newTerms.erase(foundIt);
                    newTerms.insert(newt);
                }
                continue;
            }
            // If it wasn't there, and it isn't obviously 0, insert it.
            else if (newt.factor != 0) 
                newTerms.insert(newt);
        }
    }
    p.terms = newTerms;
}

void modBy(Poly& p, Poly& q) {
    // take the first non-null term in q, write it as a sum of the others. 
}


int main() {
    using std::cout, std::endl;
    Poly p{};
    multByS(p, 12);
    multByS(p, 12);
    multByS(p, 12);
    multByS(p, 12);
    multByS(p, 12);
    cout << p.terms.size() << endl;
    
    return 0;
}
