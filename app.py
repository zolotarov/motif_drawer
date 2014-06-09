from motif_drawer_temp import *

app = Flask(__name__)
app.config.from_object('config')

@app.route('/',  methods=['POST', 'GET'])
def index():
    if request.method == "POST":
        figure_output_name = motif_drawer(server = True)
        return render_template('results.html', image_svg = 
                               '%s.svg' % figure_output_name, 
                               image_png = '%s.png' % figure_output_name)
    else:
        return render_template('index.html')


if __name__ == '__main__':
    app.run(host='0.0.0.0')