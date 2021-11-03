struct SimulationBundleForField {
  Field& E;
  Field& B;
  JField& Jmesh;
  const Grid& grid;
  const mpi::CartComm &cart;
  int timestep;
  real_t dt;
};

struct SimulationBundleForParticle {
  map<PtcArray>& particles;
  JField& J;
  std::vector<Particle<R, S>>* new_ptc_buf;
  const map<Properties>& properties;
  const Field& E;
  const Field& B;
  const Grid& grid;
  const Ensemble* ens;
  real_t dt;
  int timestep;
  util::Rng<R>& rng;
};

template <int DGrid>
struct ActionBase : public array<Range, DGrid> {
 private:
  std::string _name = "Unknown";

 public:
  ActionBase& setName(std::string name) {
    _name = name;
    return *this;
  }

  const auto& name() const noexcept { return _name; }

  virtual ~ActionBase() = default;
  // no Clone
};

template <typename Real, int DGrid, typename RealJ>
struct FieldAction : public ActionBase<DGrid> {
 public:
  virtual ~Action() = default;

  virtual void operator()(const SimulationBundleForField& bundle) const = 0;
};


template <int DGrid, typename R, template <typename> class S, typename RJ>
struct PtcAction : public ActionBase<DGrid> {
  virtual ~Action() = default;

  virtual void operator()(const SimulationBundleForParticle& bundle) const = 0;
};


struct SimulationBuilder {
 public:
  template <typename ConcreteFieldAction, typename... Args>
  ConcreteFieldAction& add_field_action(Args&&... ctor_args) {
    static_assert(std::is_base_of_v<FieldAction, ConcreteFieldAction>);
    auto& uniq_ptr = m_fld_actions.emplace_back(new ConcreteFieldAction(std::forward<Args>(ctor_args)...));
    return static_cast<ConcreteFieldAction&>(*uniq_ptr);
  }

  auto& add_particle_action();

  auto& add_extra_init( std::function<R()> ); // for RTD

  auto& set_prior_export();

  auto& add_exporter();

  auto& set_post_export();

  auto& add_custom_step();  // for checking vitals

 private:
  // these are with global ranges;
  std::vector<std::unique_ptr<FieldAction>> m_fld_actions;
  std::vector<std::unique_ptr<PtcAction>> m_ptc_actions;
};

void handle_cmd_args();

struct ConfFile {
public:
  static ConfFile load( const std::string& file );


  ConfFile operator[]( const std::string& entry ) {
    ConfFile res;
    // c++20 has std::format for compile-time format and std::vformat for
    // runtime format
    res.m_parent_entries = std::format("{}[{}]", m_parent_entries, entry);
    res.node = node[entry];
    return res;
  }

  template <typename T>
  T as() {
    // move safe_set logic here
  }

  template <typename T>
  T as_or( T val_default ) {
    // move safe_set logic here
  }

private:
  ConfFile() = default;
  std::string m_parent_entries = "".
    toml::table node;
};


struct SampleFieldAction : public FieldAction {
  
};

constexpr real_t compile_time_const = 0;

int main(int argc, char** argv) {
  handle_cmd_args();

  auto conf = ConfFile::load("some toml file");

  auto dt = conf["dt"].as<real_t>();
  auto gamma_fd = conf["pairs"]["gamma_fd"].as<real_t>();

  SimulationBuilder builder;

  {
    auto& field_act = builder.add_field_action<SampleFieldAction>();
    {
      field_act.set_name("blahblah");
    }
  }

  // auto sim = builder.build();
  // sim.start();
  // sim.print_steps(); // prints steps in order

  return 0;
}
