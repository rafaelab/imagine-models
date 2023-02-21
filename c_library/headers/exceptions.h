class NotImplementedException : public std::logic_error
{
public:
    NotImplementedException () : std::logic_error{"Function not yet implemented."} {}
};


class FieldException : public std::logic_error
{
public:
    FieldException () : std::logic_error{"Field cannot be intialized this way."} {}
};


class DivergenceException : public std::logic_error
{
public:
    DivergenceException () : std::logic_error{"The divergence of a vectorfield can only be calculated in 3 dimensions"} {}
};
